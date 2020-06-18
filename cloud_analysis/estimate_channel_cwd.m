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
OPT_DO_PLOTS = 7;

if nargin < 2
    OPTS = ...
        [ ...   
            1; ...  % compute path gain
            1; ...  % K factor
            1; ...  % delay spread
            0; ...  % compute average CIR from data
            0; ...  % compute ntap approximation
            1; ...  % write stats file    
            1; ...  % do plots
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
    rms_delay_spread_sec = nan(NN,4);     
    mean_delay_sec = nan(NN,4); 
    mean_delay_corrcoeff = NaN;
    rms_delay_spread_corrcoeff = NaN;
    cir_duration = nan(NN,1);    
    
    % Memory for calculation of average cir
    klos = 0; 
    wla = 2*wl+1;   % size of cir avg calculation buffer
    mla = floor(wla/2)+1;  %mid-point of cir avg calculation buffer
    cir_sum = zeros(wla,1);
    num_cirsum = 0;
    
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
    % LOCATION 
    %
    loc = cir_file.IQdata_CloudLocations_m_num;
    [xpos, ypos] = cloud_positions_to_records(meta, loc);
    pos = sqrt(xpos.^2 + ypos.^2);
    if isempty(pos)
        continue
    end

    % 
    % Loop through all records within the file
    %
    for kk = 1:NN
        
        % range of the CIR data record
        % r = cir_file.IQdata_Range_m(kk,3);
        r = pos(kk);
        
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
            [mean_delay_sec(kk,4), rms_delay_spread_sec(kk,4), cir_duration(kk)] = ...
                compute_delay_spread(Ts, cir);
            mean_delay_sec(kk,1:3) = [xpos(kk) ypos(kk) r];
            rms_delay_spread_sec(kk,1:3) = [xpos(kk) ypos(kk) r];
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
                
                num_cirsum = num_cirsum + 1;
                cir_sum(1:lcir0) = cir_sum(1:lcir0) + cir0;
                    
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
        pg_nans = isnan(path_gain_dB);
        r(pg_nans) = [];
        path_gain_dB(pg_nans) = [];
        
        %
        % Analyze path loss versus acquisition
        %
        if OPTS(OPT_DO_PLOTS)
        h = figure();
        yyaxis left
        plot(1:length(r(~isnan(r))),path_gain_dB(~isnan(r)))
        ylabel('Path Gain (dB)')
        yyaxis right
        plot(1:length(r(~isnan(r))),r(~isnan(r)),'.')
        ylabel('Distance from Ref. Point (m)')
        xlabel('Acquisition')
        setCommonAxisProps()    
        drawnow
        savefig(h, [fig_dir '\' mat_fname(1:end-4) '__plvacq.fig']);
        setFigureForPrinting();
        print(h,[png_dir '\' mat_fname(1:end-4) '__plvacq.png'],'-dpng','-r300')
        close(h) 
        end
        
    end % if OPTS(OPT_PATH_GAIN)
    

    if OPTS(OPT_DELAY_SPREAD) && sum(~isnan(mean_delay_sec(:,4))) > 100 && OPTS(OPT_DO_PLOTS)
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
        
        %
        % Do plot of range versus mean delay
        %
        X=mean_delay_sec(:,3); Y=mean_delay_sec(:,4);
        Z = corrcoef(X,Y); mean_delay_corrcoeff = Z(1,2);
        h = figure();
        scatter( mean_delay_sec(:,3), 1e9*mean_delay_sec(:,4) );
        xlabel('Range from Origin (m)');
        ylabel('Mean Delay (ns)');
        title('Mean Delay vs. Range from Origin')
        setCommonAxisProps();
        savefig(h,[fig_dir '\' mat_fname(1:end-4) '__duvr.fig']);
        setFigureForPrinting();
        print(h,[png_dir '\' mat_fname(1:end-4) '__duvr.png'],'-dpng','-r300')    
        close(gcf)        
        
        %
        % Do plot of range versus RMS delay
        %
        X=rms_delay_spread_sec(:,3); Y=rms_delay_spread_sec(:,4);
        Z = corrcoef(X,Y); rms_delay_spread_corrcoeff = Z(1,2);
        h = figure();
        scatter( rms_delay_spread_sec(:,3), 1e9*rms_delay_spread_sec(:,4) );
        xlabel('Range from Origin (m)');
        ylabel('RMS Delay (ns)');
        title('RMS Delay vs. Range from Origin')
        setCommonAxisProps();
        savefig(h,[fig_dir '\' mat_fname(1:end-4) '__dsvr.fig']);
        setFigureForPrinting();
        print(h,[png_dir '\' mat_fname(1:end-4) '__dsvr.png'],'-dpng','-r300')    
        close(gcf)            
        
    end %if OPTS(OPT_DELAY_SPREAD)

    if OPTS(OPT_KFACTOR) && OPTS(OPT_DO_PLOTS)
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
        h = figure(); 
        cir_avg = cir_sum/num_cirsum;

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
        savefig(h,[fig_dir '\' mat_fname(1:end-4) '__avgcir.fig']);
        setFigureForPrinting();
        print(h,[png_dir '\' mat_fname(1:end-4) '__avgcir.png'],'-dpng','-r300')    
        close(h)   

        
    end % OPTS(OPT_AVGCIR)
    
    stats = []; %#ok<NASGU>
    if OPTS(OPT_WRITE_STATS)
        % save the metrics
        stats = struct(...
            'meta',meta,...
            'path_gain_range_m', path_gain_range_m, ...
            'path_gain_dB',path_gain_dB,...
            'K',K,...
            'los', LOS, ...
            'rms_delay_spread_sec',rms_delay_spread_sec, ...
            'rms_delay_spread_corrcoeff',rms_delay_spread_corrcoeff, ...
            'mean_delay_sec',mean_delay_sec, ...
            'mean_delay_spread_corrcoeff',mean_delay_corrcoeff, ...
            'cir_duration_sec',cir_duration, ...
            'avg_cir_st', cir_avg_st); %#ok<NASGU>

%         save([stats_dir '\' mat_fname(1:end-4) '__channel_stats.mat'], 'stats')    
        save([stats_dir '\' mat_fname '__channel_stats.mat'], 'stats')    

%         % save the stats for this measurement run in local cell array
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
            nanmean(stats.K), nanmin(stats.K), nanmax(stats.K), ...
            1e9*nanmean(stats.mean_delay_sec(:,4)), 1e9*nanmin(stats.mean_delay_sec(:,4)), 1e9*nanmax(stats.mean_delay_sec(:,4)), mean_delay_corrcoeff, ...
            1e9*nanmean(stats.rms_delay_spread_sec(:,4)), 1e9*nanmin(stats.rms_delay_spread_sec(:,4)), 1e9*nanmax(stats.rms_delay_spread_sec(:,4)), rms_delay_spread_corrcoeff ...
        };
        Cstats_ii = Cstats_ii + 1;
    
    end

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
% if OPTS(OPT_PATH_GAIN)
%     cmp_pl_poly( '*_stats.mat', '.\stats', '.\figs', '.\png' )
% end

% create the delay profile files for RF emulator
if OPTS(OPT_NTAP_APPROX)
    stats2rfnestdp( '*_stats.mat', '.\stats', '.\emu' )
end

end % function

function writeStatsToFile(X)
    if isempty(X)
        return
    end
    M = cell2table(X);
    file_path = 'stats\stats.xlsx';
    VariableNames = ...
        {'Run', 'Freq', 'RX_Ant_Type', 'RX_Ant_Gain', 'TX_Ant_Type', 'TX_Ant_Gain', ...
        'Mean_K', 'Min_K', 'Max_K', ...
        'Mean_Delay_ns', 'Min_Delay_ns', 'Max_Delay_ns', 'Delay_Range_CorrCoef'...
        'Mean_Delay_Spread_ns', 'Min_Delay_Spread_ns', 'Max_Delay_Spread_ns', 'Delay_Spread_Range_CorrCoef'};
    M.Properties.VariableNames = VariableNames;
    
    M0 = [];
    if exist(file_path,'file')
        M0 = readtable(file_path);
    end
    if ~isempty(M0)
        M = [M0;M];
    end
    M.Properties.VariableNames = VariableNames;
%     writetable(M, file_path, 'Delimiter', '\t');
    writetable(M, file_path);
    
end


function b = testForStatsFile(mat_fname)
    b = false;
    if exist(['stats/' mat_fname '__channel_stats.mat'],'file')
        b = true;
    end
end


function [xpos, ypos] = cloud_positions_to_records(meta, loc)

    IQdata_CloudLocations_m_num = loc; %DATA.IQdata_CloudLocations_m_num;
    Strct_Metadata = meta; %DATA.Strct_Metadata;

    % Get x Positions and repeat each element by the number of records per acquisition
    xpos = repelem(IQdata_CloudLocations_m_num.xPositions,Strct_Metadata.NumberRecordperAcqusition_num,1);
    % Get y Positions and repeat each element by the number of records per acquisition
    ypos = repelem(IQdata_CloudLocations_m_num.yPositions,Strct_Metadata.NumberRecordperAcqusition_num,1);

end



