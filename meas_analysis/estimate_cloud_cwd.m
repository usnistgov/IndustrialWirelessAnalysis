function estimate_cloud_cwd(pattern, OPTS, TEST_DATA)
% Analyze complex impulse responses from cloud measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

TESTING = false;
if nargin == 3
    TESTING = true;
end

OPT_PATH_GAIN = 1;
OPT_DELAY_SPREAD = 2;
OPT_AVGCIR_NTAP = 3;
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
    num_los = 0;
    
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
        continue;
    end
    % plot the x,y positions
    NX = length(cir_file.IQdata_CloudLocations_m_num.xPositions);
    NY = length(cir_file.IQdata_CloudLocations_m_num.yPositions);
    h=figure();
    plot(cir_file.IQdata_CloudLocations_m_num.xPositions/lambda_m, ...
        cir_file.IQdata_CloudLocations_m_num.yPositions/lambda_m, '+')
    xlabel('xpos (\lambda)')
    ylabel('ypos (\lambda)')
    drawnow
    savefig(h,[fig_dir '\' mat_fname(1:end-4) '__position.fig']);
    setFigureForPrinting();
    print(h,[png_dir '\' mat_fname(1:end-4) '__position.png'],'-dpng','-r300')    
    close(gcf)      
    
    %
    % TODO: compute the sample variance across all of the CIR's
    %
    % show the mean and variation of h(k) at each sample time??? 
    
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
        [k_sel, nf] = select_cir_samples(cir);
        if isempty(k_sel)
            continue
        else
            if length(k_sel) < 8
                continue
            end
        end
        USE(kk) = 1;
        
        % extract the data at selected indicies
        t_k = t(k_sel);
        cir_k = cir(k_sel);

        % Compute the path loss in the cir
        % note that the CIR contains antenna gains.  We must remove the
        % bulk antenna gains using the assumption that the gain is
        % applied equally to all multi-path components.  We know that
        % this is not the true case, but without ray-tracing it is the
        % only option.
        if OPTS(OPT_PATH_GAIN)
            path_gain_dB(kk) = compute_path_gain(cir_k, ...
                TransmitterAntennaGain_dBi, ...
                ReceiverAntennaGain_dBi);           
        end
        
        % Aggregate the sums for later computation of avg CIR
        % LOS and NLOS are considered as separate classes of CIR's
        if OPTS(OPT_AVGCIR_NTAP)
            [K(kk), LOS(kk), klos] = compute_k_factor(t_k, cir_k, ns);
            if ~isnan(klos)
                inds = mlos-klos+1:wla-klos;
                num_los = num_los + 1;
                cir_sum(inds) = cir_sum(inds) + cir;
            end
        end        

        % compute delay spread parameters of the CIR 
        % because of the wrapping of energy in the FFT-based
        % correlation we must remove the trailing edge.
        if OPTS(OPT_DELAY_SPREAD)
            [mean_delay_sec(kk), rms_delay_spread_sec(kk), cir_duration(kk)] = ...
                compute_delay_spread(t_k, cir_k, nf);
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
    setCommonAxisProps()    
    xlabel('distance (\lambda)')
    ylabel('Gain (dB)')
    drawnow
    savefig(h, [fig_dir '\' mat_fname(1:end-4) '__pl.fig']);
    setFigureForPrinting();
    print(h,[png_dir '\' mat_fname(1:end-4) '__pl.png'],'-dpng','-r300')
    close(gcf)    
    
    h = figure();
    [pl_counts,pl_centers] = hist(path_gain_dB,500);
    pl_probs = cumsum(pl_counts/sum(pl_counts));
    plot(pl_centers, pl_probs, 'k');
    ylim([0 1]);
    str = 'Path Gain, G (dB)';xlabel(str,'Interpreter','Latex')
    str = 'Pr. $$ \hat{G} < G $$'; ylabel(str,'Interpreter','Latex'); 
    setCommonAxisProps();
    drawnow
    savefig(h,[fig_dir '\' mat_fname(1:end-4) '__plcdf.fig']);
    setFigureForPrinting();
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
        setCommonAxisProps();
        %title({'CDF of Average Delay', strrep(mat_fname,'_','-')})
        drawnow
        savefig(h,[fig_dir '\' mat_fname(1:end-4) '__du.fig']);
        setFigureForPrinting();
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
        setCommonAxisProps();
        drawnow
        savefig(h,[fig_dir '\' mat_fname(1:end-4) '__ds.fig']);
        h.PaperPositionMode = 'auto';
        setFigureForPrinting();
        print(h,[png_dir '\' mat_fname(1:end-4) '__ds.png'],'-dpng','-r300')    
        close(gcf)
        
    end %if OPTS(OPT_DELAY_SPREAD)
    
    % approximate an N-tap CIR from the measured CIR's
    if OPTS(OPT_AVGCIR_NTAP)
        
        NtapApprox_N = 13;
        
        cir_class_name = cir_class(cir_class_ii);
        h = figure(); 
        cir_avg = cir_sum/num_los;

        % now remove leading zeros
        cir_avg(1:find(cir_avg, 1,'first')-1) = [];
        cir_avg = cir_avg(1:wl);
        cir_avg = select_for_avg_cir(cir_avg);
        cir_avg = cir_avg/max(abs(cir_avg));
        Ncir_avg = length(cir_avg);
        t_ciravg = t(1:Ncir_avg);
        cir_avg = cir_avg(1:Ncir_avg);
        [r_t,r_h,r_ph] = reduce_taps(cir_avg,NtapApprox_N);
        r_h = r_h/max(r_h);  % normalize the approximated cir

        hold off
        %subplot(4,1,1:2)
        plot(1E9*t_ciravg, abs(cir_avg)); 
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
        setCommonAxisProps();
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
            setCommonAxisProps();
            set(gca,'OuterPosition',get(gca,'OuterPosition').*[1 1 0.95 0.95]+[0.05 0.05 0 0])
        end

        cir_avg_st(cir_class_ii).time = t_ciravg(r_t+1);
        cir_avg_st(cir_class_ii).mag = r_h;
        cir_avg_st(cir_class_ii).angle = r_ph;     

        drawnow
        savefig(h,[fig_dir '\' mat_fname(1:end-4) '__avgcir_' cell2mat(cir_class_name) '.fig']);
        setFigureForPrinting();
        print(h,[png_dir '\' mat_fname(1:end-4) '__avgcir_' cell2mat(cir_class_name) '.png'],'-dpng','-r300')    
        close(gcf)   
        
    end % OPTS(OPT_AVGCIR_NTap)
    
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
if OPT_AVGCIR_NTAP
    stats2rfnestdp( '*_stats.mat', '.\stats', '.\emu' )
end

end % function

function setCommonAxisProps()

    alw = 0.75;    % AxesLineWidth
    fsz = 10;      % Fontsize
    lw = 0.75;      % LineWidth
    msz = 3.5;       % MarkerSize
    
%    grid on
    set(gca,'XGrid','on')
    set(gca,'XMinorGrid','off')
    set(gca,'YGrid','on')
    set(gca,'YMinorGrid','off')
    set(gca,'GridAlpha',0.5)
    set(gca,'MinorGridAlpha',0.4)
    set(gca,'Fontsize',fsz)
    set(gca,'LineWidth',alw);
    set(gca,'FontName','TimesRoman')
    
    % set the line properties
    hline = get(gca,'Children');
    for h = hline(:)'
        h.LineWidth = lw;
        h.MarkerSize = msz;
    end
end

function setFigureForPrinting()
    width=3; height=3;
    %set(gcf,'InvertHardcopy','on');
    set(gcf,'PaperUnits', 'inches');
    papersize = get(gcf, 'PaperSize');
    left = (papersize(1)- width)/2;
    bottom = (papersize(2)- height)/2;
    myfiguresize = [left, bottom, width, height];
    set(gcf,'PaperPosition', myfiguresize);
end

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



function [ k, nf ] = select_cir_samples( cir )
% SELECT_CIR_SAMPLES Compute the delay spread of the input CIR
%
% Outputs:
%   k is an array of the selected indices
%   nf is the noise floor
%
% Inputs:
%   r is the range in meters from transmitter to receiver
%   cir is the real or complex valued channel impulse response
%
% Time and CIR vectors must have the same length.
%
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

% linear domain thresholds
nft = 10;           % 10 dB above noise floor
pkt = 1/31.6228;    % 15 dB from peak

% eliminate the anomalous tail components (last 8 samples)
cir = cir(1:round(length(cir)*0.5));
% cir = cir(1:end-8);

% compute the magnitude squared of the cir normailzed to the peak value
cir_mag2 = (abs(cir).^2);
cir_max = max(cir_mag2);
cir_mag2 = (abs(cir).^2)/cir_max;

% select the cir sample for use in average cir estimation
nf = mean(cir_mag2(round(length(cir)*0.8):round(length(cir)*0.9)));
k = find( cir_mag2 > nf*nft );

end