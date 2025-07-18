function estimate_channel_cwd(obj, pattern, TEST_DATA)
% Analyze complex impulse responses from measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

TESTING = false;
if nargin == 3
    TESTING = true;
end

if nargin < 1
    error 'specify a file search pattern'
end

if obj.OPTS(obj.OPT_AVGCIR) == 1
    if obj.OPTS(obj.OPT_KFACTOR) == 0
        obj.OPTS(obj.OPT_KFACTOR) = 1;
        disp('average cir estimation requires k factor estimation')
        disp('enabling k factor estimation')
    end
end

more off;

%
% query the list of measurement files
%
arr_dir = '.';          % sub directory of stored cir records
files = dir(pattern);

%
% setup data for stats text file 
%
Cstats = {};
Cstats_ii = 1;

%
% Process each file in turn
%
cir_file = [];
if TESTING
    Nfiles = 1;
else
    Nfiles = length(files);
end
for fk = 1:Nfiles
    
    % use test data or the real thing
    try 
        
        if TESTING
            disp('opening TEST_DATA');
            cir_file = TEST_DATA;
            mat_fname = [TEST_DATA.Strct_Metadata.MatFile_str '.mat'];
        else
            mat_fname = files(fk).name;             
            cir_file_path = [arr_dir '\' mat_fname]; 
            if testForStatsFile(mat_fname(1:end-4)) && obj.OPTS(obj.OPT_WRITE_STATS)
                disp(['skipping ' mat_fname])
                continue;
            end
            disp(['loading file ' mat_fname '  ...']);
            cir_file = load(cir_file_path);
        end
        
    catch me
        warning('Problem reading mat file, trying again then skipping.');
        disp(me.message)         
        warning('Skipping file...');
        continue;
    end
    
    try
        meta = cir_file.Strct_Metadata;
    catch me
        warning('problem reading meta data')
        disp(me.message);
        continue;
    end
    
    % enforce mat file str is correct, some are incorrect in the files
    meta.MatFile_str = mat_fname;
    
    % determine the site
    SITE_IND_OATS    = 0;
    SITE_IND_AAPLANT    = 1;
    SITE_IND_GBURG      = 2;
    SITE_IND_STEAM      = 3;
    if ~isempty(strfind(lower(mat_fname),'aaplant'))
        SITE_IND = SITE_IND_AAPLANT;
    elseif ~isempty(strfind(lower(mat_fname),'gburg'))
        SITE_IND = SITE_IND_GBURG;
    elseif ~isempty(strfind(lower(mat_fname),'oats'))
        SITE_IND = SITE_IND_OATS;
    else
        SITE_IND = SITE_IND_STEAM;
    end
    
    % META DATA SECTION
    %disp(meta)
    Ts = (1/meta.SampleRate_MHz_num)*1e-6;  % sample rate
    wl =  meta.CodewordLength_num;          % codeword length
    ns = meta.PNOversample_num;             % oversample rate
    t = (0:Ts:Ts*(wl-1));% -Ts*pn_over;     % time array over a burst transmission
    
    %NN = apf*rpa;
    NN = size(cir_file.IQdata,2);
    
    % truncate the last section so not to emphasize the end of the run
    NN = round(NN*0.9); 

    % Initialize memory for the metrics
    USE = nan(NN,1);  
    K = nan(NN,1);  
    LOS = nan(NN,1);  
    path_gain_dB = nan(NN,1);
    path_gain_range_m = nan(NN,1);
    path_gain_dB_poly = {};
    rms_delay_spread_sec = nan(NN,1);     
    mean_delay_sec = nan(NN,1); 
    cir_duration = nan(NN,1);    
    
    % Memory for calculation of average cir
    klos = 0; 
    wla = 2*wl+1;   % size of cir avg calculation buffer
    mla = floor(wla/2)+1;  %mid-point of cir avg calculation buffer
    cir_sum_los = zeros(wla,1);
    num_los = 0;
    cir_sum_nlos = zeros(wla,1);
    num_nlos = 0;
    
    cir_class = {'los','nlos'};
    for cir_class_ii = 1:length(cir_class)
        cir_avg_st(cir_class_ii).class = cir_class(cir_class_ii); %#ok<*AGROW>
        cir_avg_st(cir_class_ii).time = [];
        cir_avg_st(cir_class_ii).mag = [];
        cir_avg_st(cir_class_ii).angle = [];
    end

    % setup output directories
    stats_dir = 'stats';
    fig_dir = 'figs';
    png_dir = 'png';
    if ~exist(stats_dir,'dir'), mkdir(stats_dir), end;
    if ~exist(fig_dir,'dir'),   mkdir(fig_dir), end;
    if ~exist(png_dir,'dir'),   mkdir(png_dir), end; 

    %
    % ANTENNA DATA
    % 
    TransmitterAntennaGain_dBi = meta.TransmitterAntennaGain_dBi_num;
    ReceiverAntennaGain_dBi = meta.ReceiverAntennaGain_dBi_num;

    % metrics table 
    metrics_arr = [];
    reduced_cir_arr = [];
    m_kk = 1;

    % 
    % Loop through all records within the file
    %
    for kk = 1:NN
        
        % range of the CIR data record
        r = cir_file.IQdata_Range_m(kk,3);
        
        % extract the CIR for this record from the data file
        cir = cir_file.IQdata(:,kk);
        
        % compute the magnitude of the CIR samples
        cir_mag2 = abs(cir).^2;
        
        % ignore record if it is empty or the length is less than the
        % expected codeword length.  This indicates that something went
        % wrong with the instrumentation.
        if isempty(cir_mag2)
            continue;
        elseif length(cir_mag2) < wl
            continue;
        end     

        % select the sample of the cir that meet threshold criteria
        % also compute the noise floor
        [k_sel, ~, cir, pk_pwr] = obj.select_cir_samples(r, cir);
        if isempty(k_sel)
            continue
        end
        USE(kk) = 1;

        % Compute the path loss in the cir
        % note that the CIR contains antenna gains.  We must remove the
        % bulk antenna gains using the assumption that the gain is
        % applied equally to all multi-path components.  We know that
        % this is not the true case, but without ray-tracing it is the
        % only option.
        if obj.OPTS(obj.OPT_PATH_GAIN)
            path_gain_range_m(kk) = r;
            path_gain_dB(kk) = obj.compute_path_gain(cir, ...
                TransmitterAntennaGain_dBi, ...
                ReceiverAntennaGain_dBi); 
        end           

        % record the path gain
        metrics_arr(m_kk,obj.PathGain) = path_gain_dB(kk);

        % grab a reduced tap cir
        cir_mag2_red = abs(cir(1:256)).^2;
        % [cir_red_tap_t,cir_red_tap_h,cir_red_tap_ph] = reduce_taps(cir_mag2_red, obj.NtapApprox_N);
        % cir_red_tap_h = cir_mag2_red;
        cir_red_tap_h=resample(cir(1:256),1,4);
        reduced_cir_arr(m_kk,:) = cir_red_tap_h;

        % record the X and Y coordinates
        coord_x = cir_file.IQdata_Range_m(kk,4);
        coord_y = cir_file.IQdata_Range_m(kk,5);
        metrics_arr(m_kk,obj.CoordX) = coord_x;
        metrics_arr(m_kk,obj.CoordY) = coord_y;        
        
        % compute delay spread parameters of the CIR 
        % because of the wrapping of energy in the FFT-based
        % correlation we must remove the trailing edge.
        if obj.OPTS(obj.OPT_DELAY_SPREAD)
            [mean_delay_sec(kk), rms_delay_spread_sec(kk), cir_duration(kk)] = ...
                obj.compute_delay_spread(Ts, cir);
            metrics_arr(m_kk,obj.MeanDelay) = mean_delay_sec(kk);
            metrics_arr(m_kk,obj.RMSDelaySpread) = rms_delay_spread_sec(kk);
            metrics_arr(m_kk,obj.MaxDelay) = cir_duration(kk);
        end        

        % compute the K factor assuming Rician channel
        % we only consider components 10 dB above the noise floor and
        % where the peak occurs within 8 samples of beginning of the CIR 
        if obj.OPTS(obj.OPT_KFACTOR) || obj.OPTS(obj.OPT_AVGCIR)
            [K(kk), LOS(kk), k_pks] = obj.compute_k_factor(cir, ns);
            metrics_arr(m_kk,obj.RicianK) = K(kk);
            metrics_arr(m_kk,obj.LOS) = LOS(kk);
        end
        
        % Aggregate the sums for later computation of avg CIR
        % LOS and NLOS are considered as separate classes of CIR's
        % it is assumed that the clock synchronization between TX and RX is
        % working correctly
        if obj.OPTS(obj.OPT_AVGCIR)
            if pk_pwr > -100 && ~isempty(k_pks)
                cir0 = cir/max(abs(cir)); 
                if (k_pks(1)-ns > -1)
                    cir0start = k_pks(1)-ns+1;
                else
                    cir0start = 1;
                end
                cir0stop = min([k_pks(end)+ns, length(cir0)]);
                cir0 = cir0(cir0start:cir0stop);
                lcir0 = length(cir0);
                if LOS(kk) == 1 
                    num_los = num_los + 1;
                    cir_sum_los(1:lcir0) = cir_sum_los(1:lcir0) + cir0;                      
                elseif LOS(kk) == -1
                    num_nlos = num_nlos + 1;
                    cir_sum_nlos(1:lcir0) = cir_sum_nlos(1:lcir0) + cir0;                
                end      
            end   
        end      
       
        m_kk = m_kk + 1;
    end
    
    %
    % determine if the run produced enough data to form estimates.  The
    % selection of threshold was chosen arbitrarily
    %
    if sum(~isnan(USE)) < 100
        warning 'not enough CIRs passed selection to form metrics'
        disp(mat_fname)
        continue
    end

    % write the ai metrics to file
    ai_metrics_fname = [stats_dir '/' mat_fname(1:end-4) '_aimetrics.xlsx'];
    metrics_tbl_colnames = { ...
        'CoordX', 'CoordY', 'LOS', 'RicianK', ...
        'RMSDelaySpread', 'MeanDelay', 'MaxDelay', 'PathGain'};   
    metrics_tbl_vartypes = {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double' };
    metrics_tbl = array2table(metrics_arr, 'VariableNames',metrics_tbl_colnames);
    writetable(metrics_tbl, ai_metrics_fname);

    % write the reduced cirs to file
    red_cir_fname = [stats_dir '/' mat_fname(1:end-4) '_redcirs.xlsx'];
    reduced_cir_arr = [metrics_arr(:,1:8) reduced_cir_arr];
    red_cir_tbl = array2table(reduced_cir_arr);
    red_cir_tbl.Properties.VariableNames(1:8) = metrics_tbl_colnames;
    writetable(red_cir_tbl, red_cir_fname);
    
    %
    % Extract range data for off-line analysis
    %
    
    if obj.OPTS(obj.OPT_PATH_GAIN)
        
        rA = path_gain_range_m;
        r = rA;        
        pg_nans = isnan(path_gain_dB);
        r(pg_nans) = [];
        path_gain_dB(pg_nans) = [];
        
        %
        % Analyze path loss versus acquisition
        %
        if obj.OPTS(obj.OPT_DO_PLOTS)
        h = figure();
        yyaxis left
        plot(1:length(r(~isnan(r))),path_gain_dB(~isnan(r)))
        ylabel('Path Gain (dB)')
        yyaxis right
        plot(1:length(r(~isnan(r))),r(~isnan(r)),'-')
        ylabel('Distance (m)')
        xlabel('Acquisition')
        setCommonAxisProps()    
        drawnow
        savefig(h, [fig_dir '\' mat_fname(1:end-4) '__plvacq.fig']);
        setFigureForPrinting();
        print(h,[png_dir '\' mat_fname(1:end-4) '__plvacq.png'],'-dpng','-r300')
        close(h) 
        end
        
        %
        % Analyze path loss versus distance
        %
        
        R = r;
        G = path_gain_dB;
        G = G(~isnan(R));
        R = R(~isnan(R));
        
        % new way for fit
        [fr, gof] = reporting.createPwLFit(log10(R), G);
        p1 = [fr.a1, fr.b1];
        p2 = [fr.a2, fr.b2];
        %x_intersect = 10^(fr.Q);
        x_intersect = fzero(@(x) polyval(p1-p2,x),10^(fr.Q));
        path_gain_dB_poly{1} = p1;
        path_gain_dB_poly{2} = p2;
        d_val1 = log10(logspace(log10(min(R)), x_intersect, 5));
        g_val1 = polyval(p1, d_val1); 
        d_val2 = log10(logspace(x_intersect,log10(max(R)),5));
        g_val2 = polyval(p2, d_val2);   
        
        % frii as reference
        d_frii = logspace(log10(min(R)), log10(max(R)), 5);
        ff = meta.Frequency_GHz_num;
        c = physconst('LightSpeed');
        frii_fspl_dB = 10*log10(d_frii.^2) + 20*log10(ff) + 20*log10(1e9) + 20*log10(4*pi/c);         
        
        % plot the gains
        if obj.OPTS(obj.OPT_DO_PLOTS)
        h = figure();
        semilogx(R,G, 'color', [0,0,0]+0.7, 'marker', '.', 'linestyle' , 'none')
        hold on
        semilogx(10.^d_val1,g_val1,'bx-')
        semilogx(10.^d_val2,g_val2,'bo-')
        semilogx(d_frii, -frii_fspl_dB, 'k-') 
        hold off
        legend('data',sprintf('fit1(n=%.1f)',abs(p1(1)/10)),sprintf('fit2(n=%.1f)', abs(p2(1)/10)),'fspl')
        legend('Location','Best')
        
        setCommonAxisProps()    
        xlabel('distance (m)')
        ylabel('Path Gain (dB)')
        drawnow
        savefig(h, [fig_dir '\' mat_fname(1:end-4) '__pl.fig']);
        setFigureForPrinting();
        print(h,[png_dir '\' mat_fname(1:end-4) '__pl.png'],'-dpng','-r300')
        close(h) 
        end
        
               
    end % if OPTS(obj.OPT_PATH_GAIN)
    

    if obj.OPTS(obj.OPT_DELAY_SPREAD) && sum(~isnan(mean_delay_sec)) > 100 && obj.OPTS(obj.OPT_DO_PLOTS)
        %
        % Analyze the average delay of the CIR's 
        %
        X = obj.deleteoutliers(1e9*mean_delay_sec(~isnan(mean_delay_sec)));
        h = plotPDFCDF(X, 100, 'mean delay', 'ns');
        setCommonAxisProps();
        savefig(h,[fig_dir '\' mat_fname(1:end-4) '__du.fig']);
        setFigureForPrinting();
        print(h,[png_dir '\' mat_fname(1:end-4) '__du.png'],'-dpng','-r300')    
        close(gcf)

        %
        % Analyze the delay spread of the CIR's 
        %
        X = obj.deleteoutliers(1e9*rms_delay_spread_sec(~isnan(rms_delay_spread_sec)));
        h = plotPDFCDF(X, 100, 'rms delay spread', 'ns');
        setCommonAxisProps();
        savefig(h,[fig_dir '\' mat_fname(1:end-4) '__ds.fig']);
        setFigureForPrinting();
        print(h,[png_dir '\' mat_fname(1:end-4) '__ds.png'],'-dpng','-r300')    
        close(gcf)
        
    end %if OPTS(obj.OPT_DELAY_SPREAD)

    if obj.OPTS(obj.OPT_KFACTOR) && obj.OPTS(obj.OPT_DO_PLOTS)
        % 
        % Analyze the Rician K Factor Estimates CDF
        % 
        X = obj.deleteoutliers(K(~isinf(K)));
        h = plotPDFCDF(X, 300, 'Rician K Factor', 'dB');
        setCommonAxisProps();
        savefig(h,[fig_dir '\' mat_fname(1:end-4) '__Kcdf.fig']);
        setFigureForPrinting();
        print(h,[png_dir '\' mat_fname(1:end-4) '__Kcdf.png'],'-dpng','-r300')    
        close(gcf)        
    
    end % OPTS(obj.OPT_KFACTOR)
    
    % approximate an N-tap CIR from the measured CIR's
    if obj.OPTS(obj.OPT_AVGCIR) 
        
        NtapApprox_N = obj.NtapApprox_N;
        for cir_class_ii = 1:length(cir_class)
        
            cir_class_name = cir_class(cir_class_ii);
            h = figure();  
            if strcmp(cir_class_name,'los')
                cir_avg = cir_sum_los/num_los;
            else
                cir_avg = cir_sum_nlos/num_nlos;
            end
            
            % now remove leading zeros
            cir_avg(1:find(cir_avg, 1,'first')-1) = [];
            cir_avg = cir_avg(1:wl);
            cir_avg = cir_avg/max(abs(cir_avg));
            Ncir_avg = length(cir_avg);
            t_ciravg = t(1:Ncir_avg);
            cir_avg = cir_avg(1:Ncir_avg);
            [r_t,r_h,r_ph] = reduce_taps(cir_avg,NtapApprox_N);
            r_h = r_h/max(r_h);  % normalize the approximated cir

            hold off
            if obj.OPTS(obj.OPT_NTAP_APPROX)
                plot(1E9*t_ciravg, abs(cir_avg)); 
            else
                stem(1E9*t_ciravg, abs(cir_avg)); 
            end
            str = '$$\mid{h(t)}\mid$$';ylabel(str, 'Interpreter', 'Latex')
            xlabel('time (ns)')
            xlim([-20 1000]);  
            ylim([0 1.1])
            if obj.OPTS(obj.OPT_NTAP_APPROX)
                hold on; 
                stem(1E9*t_ciravg(r_t+1), abs(r_h),'d-');
                legend('Avg CIR',[num2str(NtapApprox_N) '-tap approx.'])
                hold off
            end
            reporting.setCommonAxisProps();
            set(gca,'OuterPosition',get(gca,'OuterPosition').*[1 1 0.95 0.95]+[0.05 0.05 0 0])

            cir_avg_st(cir_class_ii).time = t_ciravg(r_t+1);
            cir_avg_st(cir_class_ii).mag = r_h;
            cir_avg_st(cir_class_ii).angle = r_ph;     

            drawnow
            savefig(h,[fig_dir '\' mat_fname(1:end-4) '__avgcir_' cell2mat(cir_class_name) '.fig']);
            reporting.setFigureForPrinting(gcf);
            print(gcf,[png_dir '\' mat_fname(1:end-4) '__avgcir_' cell2mat(cir_class_name) '.png'],'-dpng','-r300')    
            close(h)   
        
        end

        
    end % OPTS(obj.OPT_AVGCIR)
    
    stats = []; %#ok<NASGU>
    if obj.OPTS(obj.OPT_WRITE_STATS)
        % save the metrics
        stats = struct(...
            'meta',meta,...
            'path_gain_range_m', path_gain_range_m, ...
            'path_gain_dB',path_gain_dB,...
            'path_gain_dB_poly',path_gain_dB_poly{2},...
            'K',K,...
            'los', LOS, ...
            'rms_delay_spread_sec',rms_delay_spread_sec, ...
            'mean_delay_sec',mean_delay_sec, ...
            'cir_duration_sec',cir_duration, ...
            'avg_cir_st', cir_avg_st); %#ok<NASGU>

        save([stats_dir '\' mat_fname(1:end-4) '__channel_stats.mat'], 'stats')    
    
    end

    % explicit clear of large memory
    cir_file = [];  %#ok<NASGU>

end

% add entry to the stats text file
if obj.OPTS(obj.OPT_WRITE_STATS)
    writeStatsToFile(Cstats);
end

%
% Create summary data
%

% create the aggregate polynomial for path loss
if obj.OPTS(obj.OPT_PATH_GAIN)
    obj.cmp_pl_poly( '*_stats.mat', '.\stats', '.\figs', '.\png' )
end

% create the delay profile files for RF emulator
if obj.OPTS(obj.OPT_NTAP_APPROX)
    obj.stats2rfnestdp( '*_stats.mat', '.\stats', '..\emu' )
end

end % function

function writeStatsToFile(X)
    if isempty(X)
        return
    end
    M = cell2table(X);
    file_path = '../stats.dat';
    VariableNames = ...
        {'Run', 'Freq', 'RX_Ant_Type', 'RX_Ant_Gain', 'TX_Ant_Type', 'TX_Ant_Gain', ...
        'Path_Gain_Poly_Slope', 'Path_Gain_Poly_YInt', 'Mean_K', 'Min_K', 'Max_K', ...
        'Mean_Delay_ns', 'Min_Delay_ns', 'Max_Delay_ns', ...
        'Mean_Delay_Spread_ns', 'Min_Delay_Spread_ns', 'Max_Delay_Spread_ns'};
    M.Properties.VariableNames = VariableNames;
    
    M0 = [];
    if exist(file_path,'file')
        M0 = readtable(file_path);
    end
    if ~isempty(M0)
        M = [M0;M];
    end
    M.Properties.VariableNames = VariableNames;
    writetable(M, file_path, 'Delimiter', '\t');
    
end


function b = testForStatsFile(mat_fname)
    b = false;
    if exist(['stats/' mat_fname '__channel_stats.mat'],'file')
        b = true;
    end
end






