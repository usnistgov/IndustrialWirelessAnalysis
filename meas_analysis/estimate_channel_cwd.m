function estimate_channel_cwd(pattern, OPTS, TEST_DATA)
% Analyze complex impulse responses from measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

TESTING = false;
if nargin == 3
    TESTING = true;
end

OPT_PATH_GAIN = 1;
OPT_KFACTOR = 2;
OPT_DELAY_SPREAD = 3;
OPT_AVGCIR = 4;
OPT_NTAP_APPROX = 5;
OPT_WRITE_STATS = 6;

if nargin < 2
    OPTS = ...
        [ ...   
            1; ...  % compute path gain
            1; ...  % K factor
            1; ...  % delay spread
            1; ...  % compute average CIR from data
            1; ...  % compute ntap approximation
            1; ...  % write stats file    
        ]; 
end

if nargin < 1
    error 'specify a file search pattern'
end

if OPTS(OPT_AVGCIR) == 1
    if OPTS(OPT_KFACTOR) == 0
        OPTS(OPT_KFACTOR) = 1;
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
        
        mat_fname = files(fk).name; 
        
        cir_file_path = [arr_dir '\' mat_fname];       
        if TESTING
            disp('opening TEST_DATA');
            cir_file = TEST_DATA;
            mat_fname = TEST_DATA.Strct_Metadata.MatFile_str;
        else
            if testForStatsFile(mat_fname(1:end-4)) && OPTS(OPT_WRITE_STATS)
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
        [k_sel, ~, cir, pk_pwr] = select_cir_samples(r, cir);
        if isempty(k_sel)
            continue
        end
        USE(kk) = 1;
        
        % compute delay spread parameters of the CIR 
        % because of the wrapping of energy in the FFT-based
        % correlation we must remove the trailing edge.
        if OPTS(OPT_DELAY_SPREAD)
            [mean_delay_sec(kk), rms_delay_spread_sec(kk), cir_duration(kk)] = ...
                compute_delay_spread(Ts, cir);
        end        

        % Compute the path loss in the cir
        % note that the CIR contains antenna gains.  We must remove the
        % bulk antenna gains using the assumption that the gain is
        % applied equally to all multi-path components.  We know that
        % this is not the true case, but without ray-tracing it is the
        % only option.
        if OPTS(OPT_PATH_GAIN)
            path_gain_range_m(kk) = r;
            path_gain_dB(kk) = compute_path_gain(cir, ...
                TransmitterAntennaGain_dBi, ...
                ReceiverAntennaGain_dBi);           
        end

        % compute the K factor assuming Rician channel
        % we only consider components 10 dB above the noise floor and
        % where the peak occurs within 8 samples of beginning of the CIR 
        if OPTS(OPT_KFACTOR) || OPTS(OPT_AVGCIR)
            [K(kk), LOS(kk), k_pks] = compute_k_factor(cir, ns);
        end
        
        % Aggregate the sums for later computation of avg CIR
        % LOS and NLOS are considered as separate classes of CIR's
        % it is assumed that the clock synchronization between TX and RX is
        % working correctly
        if OPTS(OPT_AVGCIR)
            if pk_pwr > -100
                cir0 = cir/max(abs(cir)); 
                if (k_pks(1)-ns > -1)
                    cir0start = k_pks(1)-ns+1;
                else
                    cir0start = 1;
                end
                cir0stop = k_pks(end)+ns;
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
    
    %
    % Extract range data for off-line analysis
    %
    
    if OPTS(OPT_PATH_GAIN)
        
        rA = path_gain_range_m;
        r = rA;        
        
        %
        % Analyze path loss versus acquisition
        %
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
        
        
        %
        % Analyze path loss versus distance
        %
        
        infpt = 0;
        if SITE_IND == SITE_IND_AAPLANT
            infpt = 30;
        elseif SITE_IND == SITE_IND_GBURG
        	infpt = 10;
        elseif SITE_IND == SITE_IND_STEAM
        	infpt = 5;
        elseif SITE_IND == SITE_IND_OATS;
            infpt = 5;
        end
        R = r;
        G = path_gain_dB;
        G = G(~isnan(R));
        R = R(~isnan(R));
        
        % lower segment
        Gfit1 = G(R<infpt);
        Rfit1 = R(R<infpt);
        if ~isempty(Rfit1)
            p1 = polyfit(log10(Rfit1), Gfit1, 1);
        else
            p1 = [0 infpt];
        end
        path_gain_dB_poly{1} = p1;

        % upper segment
        Gfit2 = G(R>infpt);
        Rfit2 = R(R>infpt);
        p2 = polyfit(log10(Rfit2), Gfit2, 1);   
        path_gain_dB_poly{2} = p2;

        % find draw points
        if ~isempty(Rfit1)
            x_intersect = fzero(@(x) polyval(p1-p2,x),log10(infpt));
            %x_intersect = fzero(@(x) polyval(p1-p2,x),[log10(min(Rfit1)), log10(max(Rfit2))]);
            d_val1 = log10(logspace(log10(min(Rfit1)), x_intersect, 5));
            g_val1 = polyval(p1, d_val1);
        else
            x_intersect = log10(min(R));
        end
        d_val2 = log10(logspace(x_intersect, log10(max(Rfit2)), 5));
        g_val2 = polyval(p2, d_val2);
        
        % frii as reference
        d_frii = logspace(log10(min(R)), log10(max(R)), 5);
        ff = meta.Frequency_GHz_num;
        c = physconst('LightSpeed');
        frii_fspl_dB = 10*log10(d_frii.^2) + 20*log10(ff) + 20*log10(1e9) + 20*log10(4*pi/c);         
        
        % plot the gains
        h = figure();
        semilogx(R,G, 'color', [0,0,0]+0.7, 'marker', '.', 'linestyle' , 'none')
        hold on
        if ~isempty(Rfit1)
            semilogx(10.^d_val1,g_val1,'bx-')
        end
        semilogx(10.^d_val2,g_val2,'bo-')
        semilogx(d_frii, -frii_fspl_dB, 'k-') 
        hold off
        if ~isempty(Rfit1) && ~isempty(Rfit2)
            legend('data',sprintf('fit1(n=%.1f)',abs(p1(1)/10)),sprintf('fit2(n=%.1f)', abs(p2(1)/10)),'fspl')
        else
            legend('data',sprintf('fit(n=%.1f)',abs(p2(1)/10)),'fspl')
        end
        legend('Location','Best')
        
        setCommonAxisProps()    
        xlabel('distance (m)')
        ylabel('Path Gain (dB)')
        drawnow
        savefig(h, [fig_dir '\' mat_fname(1:end-4) '__pl.fig']);
        setFigureForPrinting();
        print(h,[png_dir '\' mat_fname(1:end-4) '__pl.png'],'-dpng','-r300')
        close(h) 
        
    end % if OPTS(OPT_PATH_GAIN)
    

    if OPTS(OPT_DELAY_SPREAD) && sum(~isnan(mean_delay_sec)) > 100
        %
        % Analyze the average delay of the CIR's 
        %
        X = deleteoutliers(1e9*mean_delay_sec(~isnan(mean_delay_sec)));
        h = plotPDFCDF(X, 100, 'mean delay', 'ns');
        setCommonAxisProps();
        savefig(h,[fig_dir '\' mat_fname(1:end-4) '__du.fig']);
        setFigureForPrinting();
        print(h,[png_dir '\' mat_fname(1:end-4) '__du.png'],'-dpng','-r300')    
        close(gcf)

        %
        % Analyze the delay spread of the CIR's 
        %
        X = deleteoutliers(1e9*rms_delay_spread_sec(~isnan(rms_delay_spread_sec)));
        h = plotPDFCDF(X, 100, 'rms delay spread', 'ns');
        setCommonAxisProps();
        savefig(h,[fig_dir '\' mat_fname(1:end-4) '__ds.fig']);
        setFigureForPrinting();
        print(h,[png_dir '\' mat_fname(1:end-4) '__ds.png'],'-dpng','-r300')    
        close(gcf)
        
    end %if OPTS(OPT_DELAY_SPREAD)

    if OPTS(OPT_KFACTOR)
        % 
        % Analyze the Rician K Factor Estimates CDF
        % 
        X = deleteoutliers(K(~isinf(K)));
        h = plotPDFCDF(X, 300, 'Rician K Factor', 'dB');
        setCommonAxisProps();
        savefig(h,[fig_dir '\' mat_fname(1:end-4) '__Kcdf.fig']);
        setFigureForPrinting();
        print(h,[png_dir '\' mat_fname(1:end-4) '__Kcdf.png'],'-dpng','-r300')    
        close(gcf)        
    
    end % OPTS(OPT_KFACTOR)
    
    % approximate an N-tap CIR from the measured CIR's
    if OPTS(OPT_AVGCIR)
        
        NtapApprox_N = 13;
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
            if OPTS(OPT_NTAP_APPROX)
                plot(1E9*t_ciravg, abs(cir_avg)); 
            else
                stem(1E9*t_ciravg, abs(cir_avg)); 
            end
            str = '$$\mid{h(t)}\mid$$';ylabel(str, 'Interpreter', 'Latex')
            xlabel('time (ns)')
            xlim([-20 1000]);  
            ylim([0 1.1])
            if OPTS(OPT_NTAP_APPROX)
                hold on; 
                stem(1E9*t_ciravg(r_t+1), abs(r_h),'d-');
                legend('Avg CIR',[num2str(NtapApprox_N) '-tap approx.'])
                hold off
            end
            setCommonAxisProps();
            set(gca,'OuterPosition',get(gca,'OuterPosition').*[1 1 0.95 0.95]+[0.05 0.05 0 0])

            cir_avg_st(cir_class_ii).time = t_ciravg(r_t+1);
            cir_avg_st(cir_class_ii).mag = r_h;
            cir_avg_st(cir_class_ii).angle = r_ph;     

            drawnow
            savefig(h,[fig_dir '\' mat_fname(1:end-4) '__avgcir_' cell2mat(cir_class_name) '.fig']);
            setFigureForPrinting();
            print(h,[png_dir '\' mat_fname(1:end-4) '__avgcir_' cell2mat(cir_class_name) '.png'],'-dpng','-r300')    
            close(h)   
        
        end

        
    end % OPTS(OPT_AVGCIR)
    
    if OPTS(OPT_WRITE_STATS)
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
            'avg_cir_st', cir_avg_st);

        save([stats_dir '\' mat_fname(1:end-4) '__channel_stats.mat'], 'stats')
    end

    % save the stats for this measurement run
    RxPol = 'U';
    if strfind(meta.ReceiverAntenna_str,'V Pol')
        RxPol = 'V';
    elseif ~isempty(strfind(meta.ReceiverAntenna_str,'Cross Pol')) || ...
            ~isempty(strfind(meta.ReceiverAntenna_str,'X Pol')) || ...
            ~isempty(strfind(meta.ReceiverAntenna_str,'Long Pol'))
        RxPol = 'X';
    end
    TxPol = 'U';
    if ~isempty(strfind(meta.TransmitterAntenna_str,'V Pol'))
        TxPol = 'V';
    elseif ~isempty(strfind(meta.TransmitterAntenna_str,'Cross Pol')) || ...
            ~isempty(strfind(meta.TransmitterAntenna_str,'X Pol')) || ...
            ~isempty(strfind(meta.TransmitterAntenna_str,'Long Pol'))
        TxPol = 'X';
    end    
    Cstats(Cstats_ii,:) = {   
        meta.MatFile_str, meta.Frequency_GHz_num, ...
        RxPol, meta.ReceiverAntennaGain_dBi_num, ...
        TxPol, meta.TransmitterAntennaGain_dBi_num, ...
        stats.path_gain_dB_poly(1), stats.path_gain_dB_poly(2), ...
        nanmean(stats.K), nanmin(stats.K), nanmax(stats.K) ...
        1e9*nanmean(stats.mean_delay_sec), 1e9*nanmin(stats.mean_delay_sec), 1e9*nanmax(stats.mean_delay_sec)...
        1e9*nanmean(stats.rms_delay_spread_sec), 1e9*nanmin(stats.rms_delay_spread_sec), 1e9*nanmax(stats.rms_delay_spread_sec) ...
    };
    Cstats_ii = Cstats_ii + 1;

    % explicit clear of large memory
    cir_file = [];  %#ok<NASGU>

end

% add entry to the stats text file
if OPTS(OPT_WRITE_STATS)
    writeStatsToFile(Cstats);
end

%
% Create summary data
%

% create the aggregate polynomial for path loss
if OPTS(OPT_PATH_GAIN)
    cmp_pl_poly( '*_stats.mat', '.\stats', '.\figs', '.\png' )
end

% create the delay profile files for RF emulator
if OPTS(OPT_NTAP_APPROX)
    stats2rfnestdp( '*_stats.mat', '.\stats', '..\emu' )
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






