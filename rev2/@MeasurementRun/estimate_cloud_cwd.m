function estimate_cloud_cwd(obj, pattern, TEST_DATA)
% Analyze complex impulse responses from cloud measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

TESTING = false;
if nargin == 3
    TESTING = true;
end

OPTS = obj.OPTS;
OPT_PATH_GAIN = 1;
OPT_DELAY_SPREAD = 2;
OPT_AVGCIR = 3;
OPT_NTAP_APPROX = 4;
OPT_WRITE_STATS = 5;
if nargin < 2
    OPTS = ...
        [ ...   
            1; ...  % compute path gain
            1; ...  % delay spread
            1; ...  % compute average CIR from data
            0; ...  % n tap approximation all cir records
            1; ...  % append stats file
        ]; 
end

if nargin < 1
    error 'specify a file search pattern'
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
        else
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
    
    % META DATA SECTION
    disp(meta)
    Ts = (1/meta.SampleRate_MHz_num)*1e-6;  % sample rate
    wl =  meta.CodewordLength_num;          % codeword length
    ns = meta.PNOversample_num;             % oversample rate
    t = (0:Ts:Ts*(wl-1));% -Ts*pn_over;     % time array over a burst transmission
    f = meta.Frequency_GHz_num;
    lambda_m = 299792458/(f*1e9);
    
    %NN = apf*rpa;
    NN = size(cir_file.IQdata,2);

    % Initialize memory for the metrics
    USE = nan(NN,1);  
    K = nan(NN,1);  
    LOS = nan(NN,1);  
    path_gain_dB = nan(NN,1);  
    rms_delay_spread_sec = nan(NN,1);     
    mean_delay_sec = nan(NN,1); 
    cir_duration = nan(NN,1);    
    
    % Memory for calculation of average cir
    wla = 2*wl+1;   % size of cir avg calculation buffer
    mlos = ceil(wla/2);  %mid-point of cir avg calculation buffer
    cir_sum = zeros(wla,1);
    cir_sum_los = zeros(wla,1);
    num_los = 0;
    cir_sum_nlos = zeros(wla,1);
    num_nlos = 0;
    
    cir_class = {'los','nlos'};
    for cir_class_ii = 1:length(cir_class)
        cir_avg_st(cir_class_ii).class = cir_class; %#ok<*AGROW>
        cir_avg_st(cir_class_ii).time = [];
        cir_avg_st(cir_class_ii).mag = [];
        cir_avg_st(cir_class_ii).angle = [];
    end

    % setup output directories
    stats_dir = 'stats';
    fig_dir = 'figs';
    png_dir = 'png';
    mkdir(stats_dir);
    mkdir(fig_dir); 
    mkdir(png_dir); 

    %
    % ANTENNA DATA
    % 
    TransmitterAntennaGain_dBi = meta.TransmitterAntennaGain_dBi_num;
    ReceiverAntennaGain_dBi = meta.ReceiverAntennaGain_dBi_num;
    
    % verify that x,y position data is available for this file
    if isempty(cir_file.IQdata_CloudLocations_m_num.xPositions)
        warning 'no position data found'
        return;
    end
    % plot the x,y positions
    CoordX = cir_file.IQdata_CloudLocations_m_num.xPositions;
    CoordY = cir_file.IQdata_CloudLocations_m_num.yPositions;
    NX = length(CoordX);
    NY = length(CoordY);
    h=figure();
    plot(CoordX/lambda_m, CoordY/lambda_m, '+')
    xlabel('xpos (\lambda)')
    ylabel('ypos (\lambda)')
    drawnow
    savefig(h,[fig_dir '\' mat_fname(1:end-4) '__position.fig']);
    reporting.setFigureForPrinting(gcf);
    print(h,[png_dir '\' mat_fname(1:end-4) '__position.png'],'-dpng','-r300')    
    close(gcf)      
    
    % metrics table 
    metrics_arr = [];
    reduced_cir_arr = [];
    m_kk = 1;
    
    % 
    % Loop through all records within the file
    %
    for kk = 1:NN
        
        % position of the CIR data record
        k_pos = mod(kk-1,NX)+1;
        POS = [ cir_file.IQdata_CloudLocations_m_num.xPositions(k_pos), ...
                cir_file.IQdata_CloudLocations_m_num.yPositions(k_pos) ];
            
        % distance from origin in wavelengths
        r(kk) = sqrt(sum(POS.^2))/lambda_m;
        
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

        % record frequency
        metrics_arr(m_kk, MeasurementRunMetric.Freq) = f;

        % record the antenna position
        % record the X and Y coordinates
        % metrics_arr(m_kk,MeasurementRunMetric.CoordX) = meta.Rx_xyz_m_cll{1};
        % metrics_arr(m_kk,MeasurementRunMetric.CoordY) = meta.Rx_xyz_m_cll{2};            
        metrics_arr(m_kk,MeasurementRunMetric.CoordX) = POS(1);
        metrics_arr(m_kk,MeasurementRunMetric.CoordY) = POS(2);

        % record the path gain
        metrics_arr(m_kk,MeasurementRunMetric.PathGain) = nan;        

        % grab a reduced tap cir
        cir_mag2_red = abs(cir(1:256)).^2;
        % [cir_red_tap_t,cir_red_tap_h,cir_red_tap_ph] = ...
        %       reduce_taps(cir_mag2_red, obj.NtapApprox_N);
        % cir_red_tap_h = cir_mag2_red;
        cir_red_tap_h=resample(cir(1:256),1,4);
        reduced_cir_arr(m_kk,:) = cir_red_tap_h;      
        
        % compute delay spread parameters of the CIR 
        % because of the wrapping of energy in the FFT-based
        % correlation we must remove the trailing edge.
        if OPTS(OPT_DELAY_SPREAD)
            if pk_pwr > -100
                [mean_delay_sec(kk), rms_delay_spread_sec(kk), cir_duration(kk)] = ...
                    obj.compute_delay_spread(Ts, cir);
                metrics_arr(m_kk,MeasurementRunMetric.MeanDelay) = mean_delay_sec(kk);
                metrics_arr(m_kk,MeasurementRunMetric.RMSDelaySpread) = rms_delay_spread_sec(kk);
                metrics_arr(m_kk,MeasurementRunMetric.MaxDelay) = cir_duration(kk);                
            end
        end    

        % compute the K factor assuming Rician channel
        % we only consider components 10 dB above the noise floor and
        % where the peak occurs within 8 samples of beginning of the CIR 
        if obj.OPTS(obj.OPT_KFACTOR) || obj.OPTS(obj.OPT_AVGCIR)
            [K(kk), LOS(kk), k_pks] = obj.compute_k_factor(cir, ns);
            metrics_arr(m_kk,MeasurementRunMetric.RicianK) = K(kk);
            metrics_arr(m_kk,MeasurementRunMetric.LOS) = LOS(kk);
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

        % Compute the path loss in the cir
        % note that the CIR contains antenna gains.  We must remove the
        % bulk antenna gains using the assumption that the gain is
        % applied equally to all multi-path components.  We know that
        % this is not the true case, but without ray-tracing it is the
        % only option.
        if OPTS(OPT_PATH_GAIN)
            path_gain_dB(kk) = obj.compute_path_gain(cir, ...
                TransmitterAntennaGain_dBi, ...
                ReceiverAntennaGain_dBi);   
            metrics_arr(m_kk, MeasurementRunMetric.PathGain) = path_gain_dB(kk);
        end
        
        % Aggregate the sums for later computation of avg CIR
        % LOS and NLOS are considered as separate classes of CIR's
        if OPTS(OPT_AVGCIR)
            [~, ~, k_pks] = obj.compute_k_factor(cir, ns);
            if pk_pwr > -100  && ~isempty(k_pks)
                cir0 = cir/max(abs(cir));
                num_los = num_los + 1;
                inds = k_pks-k_pks(1)+1;
                cir_sum(inds) = cir_sum(inds) + cir0(k_pks);
                %stem(abs(cir_sum)), xlim([0 60]), drawnow
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
    metrics_tbl = array2table(metrics_arr, 'VariableNames',obj.metrics_tbl_colnames);    
    % ai_metrics_fname = [stats_dir '/' mat_fname(1:end-4) '_aimetrics.csv'];
    % writetable(metrics_tbl, ai_metrics_fname, 'WriteMode','overwrite');
    ai_metrics_fname = [stats_dir '/' obj.AI_METRICS_FNAME_OUT];
    writetable(metrics_tbl, ai_metrics_fname, 'WriteMode','append');

    % write the reduced cirs to file
    reduced_cir_arr = [metrics_arr(:,1:obj.NBaseColumnNames) reduced_cir_arr];
    red_cir_tbl = array2table(reduced_cir_arr);
    red_cir_tbl.Properties.VariableNames(1:obj.NBaseColumnNames) = obj.metrics_tbl_colnames(1:obj.NBaseColumnNames);
    % red_cir_fname = [stats_dir '/' mat_fname(1:end-4) '_redcirs.csv'];
    % writetable(red_cir_tbl, red_cir_fname, 'WriteMode','overwrite');  
    red_cir_fname = [stats_dir '/' obj.AI_REDTAP_FNAME_OUT];
    writetable(red_cir_tbl, red_cir_fname, 'WriteMode','append');  
    
    if OPTS(OPT_PATH_GAIN)
    %
    % Analyze path loss versus distance
    %
    pl_p = path_gain_dB;
    r_p = r(~isnan(pl_p));
    pl_p = pl_p(~isnan(pl_p));
    r_p_fit = r_p(:);
    pl_p_fit = pl_p(:);
    path_gain_dB_poly = polyfit(r_p_fit,pl_p_fit,1);
    
    h = figure();
    stdPathGain = std(path_gain_dB(~isnan(path_gain_dB)));
    r_p_plot = linspace(min(r_p_fit), max(r_p_fit), 10);
    pl_poly_vals = polyval(path_gain_dB_poly, r_p_plot);
    plot(r_p, pl_p, 'color', [0,0,0]+0.4, 'marker', '.', 'linestyle' , 'none'); 
    hold on
    plot(r_p_plot, pl_poly_vals, 'k-', ...
       r_p_plot, repmat(pl_poly_vals(:),1,2)+2*stdPathGain*[ones(10,1) -ones(10,1)], ...
        'k--', 'LineWidth', 1.0);
    hold off
    legend({'measured', ...
        sprintf('%0.2fx + %0.1f',path_gain_dB_poly), ...
        '+/- 2\sigma'}, 'Location', 'best');
    reporting.setCommonAxisProps()    
    xlabel('distance (\lambda)')
    ylabel('Gain (dB)')
    drawnow
    savefig(h, [fig_dir '\' mat_fname(1:end-4) '__pl.fig']);
    reporting.setFigureForPrinting(gcf);
    print(h,[png_dir '\' mat_fname(1:end-4) '__pl.png'],'-dpng','-r300')
    close(gcf)    
    
    h = figure();
    [pl_counts,pl_centers] = hist(path_gain_dB,500);
    pl_probs = cumsum(pl_counts/sum(pl_counts));
    plot(pl_centers, pl_probs, 'k');
    ylim([0 1]);
    str = 'Path Gain, G (dB)';xlabel(str,'Interpreter','Latex')
    str = 'Pr. $$ \hat{G} < G $$'; ylabel(str,'Interpreter','Latex'); 
    reporting.setCommonAxisProps();
    drawnow
    savefig(h,[fig_dir '\' mat_fname(1:end-4) '__plcdf.fig']);
    reporting.setFigureForPrinting(gcf);
    print(h,[png_dir '\' mat_fname(1:end-4) '__plcdf.png'],'-dpng','-r300')    

    close(h)    
    
    end % if OPTS(OPT_PATH_GAIN)
    

    if OPTS(OPT_DELAY_SPREAD) && sum(~isnan(mean_delay_sec)) > 100
        %
        % Analyze the average delay of the CIR's 
        %
        h = figure();      
        [du_counts,du_centers] = hist(1e9*mean_delay_sec,50000);
        du_probs = cumsum(du_counts/sum(du_counts));
        plot(du_centers, du_probs, 'k');
        ylim([0 1]);
        xlim([0 max(du_centers(du_probs<0.995))]);
        str = 'average delay, $$\tau_D$$ (ns)';xlabel(str,'Interpreter','Latex')
        str = 'Pr. $$ \hat{\tau_D} < {\tau_D} $$'; ylabel(str,'Interpreter','Latex'); 
        reporting.setCommonAxisProps();
        %title({'CDF of Average Delay', strrep(mat_fname,'_','-')})
        drawnow
        savefig(h,[fig_dir '\' mat_fname(1:end-4) '__du.fig']);
        reporting.setFigureForPrinting(gcf);
        print(h,[png_dir '\' mat_fname(1:end-4) '__du.png'],'-dpng','-r300')    
        close(gcf)

        %
        % Analyze the delay spread of the CIR's 
        %
        h = figure();      
        [ds_counts,ds_centers] = hist(1e9*rms_delay_spread_sec,5000);
        ds_probs = cumsum(ds_counts/sum(ds_counts));
        plot(ds_centers, ds_probs, 'k');
        ylim([0 1]);
        str = 'rms delay spread, $$S$$ (ns)';xlabel(str,'Interpreter','Latex')
        str = 'Pr. $$\hat{S} < S$$'; ylabel(str,'Interpreter','Latex');
        reporting.setCommonAxisProps();
        drawnow
        savefig(h,[fig_dir '\' mat_fname(1:end-4) '__ds.fig']);
        h.PaperPositionMode = 'auto';
        reporting.setFigureForPrinting(gcf);
        print(h,[png_dir '\' mat_fname(1:end-4) '__ds.png'],'-dpng','-r300')    
        close(gcf)
        
    end %if OPTS(OPT_DELAY_SPREAD)
    
    % approximate an N-tap CIR from the measured CIR's
    if OPTS(OPT_AVGCIR)
        
        NtapApprox_N = 13;
        
        cir_class_name = cir_class(cir_class_ii);
        h = figure(); 
        cir_avg = cir_sum/num_los;

        % now remove leading zeros
        cir_avg(1:find(cir_avg, 1,'first')-1) = [];
        cir_avg = cir_avg(1:wl);
        cir_avg = obj.select_for_avg_cir(cir_avg);
        cir_avg = cir_avg/max(abs(cir_avg));
        Ncir_avg = length(cir_avg);
        t_ciravg = t(1:Ncir_avg);
        cir_avg = cir_avg(1:Ncir_avg);
        [r_t,r_h,r_ph] = reduce_taps(cir_avg,NtapApprox_N);
        r_h = r_h/max(r_h);  % normalize the approximated cir

        hold off
        %subplot(4,1,1:2)
        if OPTS(OPT_NTAP_APPROX)
            plot(1E9*t_ciravg, abs(cir_avg)); 
        else
            stem(1E9*t_ciravg, abs(cir_avg)); 
        end
        str = '$$\mid{h(t)}\mid$$';ylabel(str, 'Interpreter', 'Latex')
        %set(gca,'XTickLabel','')
        xlabel('time (ns)')
        xlim([0 1000]);         
        if OPTS(OPT_NTAP_APPROX)
            hold on; 
            stem(1E9*t_ciravg(r_t+1), abs(r_h),'d-');
            legend('Avg CIR','N-tap approx.')
            hold off
        end        
        reporting.setCommonAxisProps();
        set(gca,'OuterPosition',get(gca,'OuterPosition').*[1 1 0.95 0.95]+[0.05 0.05 0 0])

        if 0
            subplot(4,1,3:4)
            plot(1E9*t_ciravg, angle(cir_avg));
            str = '$$\angle{h(t)}$$';ylabel(str, 'Interpreter', 'Latex')
            set(gca,...
                 'ylim',[-2*pi() 2*pi()],...
                 'ytick',[-2*pi() 0 2*pi()],...
                 'yticklabel',{'-2\pi' '0' '2\pi'})
            xlim([0 1000]); 
            hold on
            stem(1E9*t_ciravg(r_t+1), r_ph, 'r')
            hold off
            reporting.setCommonAxisProps();
            set(gca,'OuterPosition',get(gca,'OuterPosition').*[1 1 0.95 0.95]+[0.05 0.05 0 0])
        end

        cir_avg_st(cir_class_ii).time = t_ciravg(r_t+1);
        cir_avg_st(cir_class_ii).mag = r_h;
        cir_avg_st(cir_class_ii).angle = r_ph;     

        drawnow
        savefig(h,[fig_dir '\' mat_fname(1:end-4) '__avgcir_' cell2mat(cir_class_name) '.fig']);
        reporting.setFigureForPrinting(gcf);
        print(h,[png_dir '\' mat_fname(1:end-4) '__avgcir_' cell2mat(cir_class_name) '.png'],'-dpng','-r300')    
        close(gcf)   
        
    end % OPTS(OPT_AVGCIR)
    
    % save the metrics
    stats = struct(...
        'meta',meta,...
        'path_gain_dB',path_gain_dB,...
        'rms_delay_spread_sec',rms_delay_spread_sec, ...
        'mean_delay_sec',mean_delay_sec, ...
        'cir_duration_sec',cir_duration, ...
        'avg_cir_st', cir_avg_st);
    
    save([stats_dir '\' mat_fname(1:end-4) '__channel_stats.mat'], 'stats')

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
        1e9*nanmean(stats.mean_delay_sec), 1e9*nanmin(stats.mean_delay_sec), 1e9*nanmax(stats.mean_delay_sec)...
        1e9*nanmean(stats.rms_delay_spread_sec), 1e9*nanmin(stats.rms_delay_spread_sec), 1e9*nanmax(stats.rms_delay_spread_sec) ...
    };
    Cstats_ii = Cstats_ii + 1;

    % explicit clear of large memory
    cir_file = [];  %#ok<NASGU>

end

% add entry to the stats text file
if OPT_WRITE_STATS
    writeStatsToFile(Cstats);
end

% create the delay profile files for RF emulator
if OPT_AVGCIR
    obj.stats2rfnestdp( '*_stats.mat', '.\stats', '.\emu' )
end

close all;

end % function

function writeStatsToFile(X)
    if isempty(X)
        return
    end
    M = cell2table(X);
    file_path = '..\cloud_stats.dat';
    VariableNames = ...
        {'Run', 'Freq', 'RX_Ant_Type', 'RX_Ant_Gain', 'TX_Ant_Type', 'TX_Ant_Gain', ...
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
