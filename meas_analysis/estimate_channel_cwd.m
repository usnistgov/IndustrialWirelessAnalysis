function estimate_channel_cwd(pattern, doall, figvis, TEST_DATA)
% Analyze complex impulse responses from measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

TESTING = false;
if nargin == 4
    TESTING = true;
end

if nargin < 3
    figvis = true;
end
if nargin < 2
    doall = false;
end

more off; 
grid minor;

%
peaks = [];
peaks_t = [];
peaks_k = [];
K = [];
cir_file = [];
arr_dir = '.';

files = dir(pattern);
for fk = 1:length(files)
    
    % check for semphore
    mat_fname = files(fk).name;
    sem_fname = ['semv\' mat_fname(1:end-4) '__sem.mat'];
    if ~doall && exist(sem_fname, 'file')
        disp('file exists. continuing...');
        continue;
    end
    try 
        cir_file_path = [arr_dir '\' mat_fname];
        disp(['opening ' mat_fname]);
        if TESTING
            cir_file = TEST_DATA;
        else
            cir_file = load(cir_file_path);
        end
    catch me
        warning('Problem reading mat file, trying again then skipping.');
        disp(me.message)
        try
        cir_file = load(mat_fname);
        catch
            disp(me.message)            
            warning('Skipping file...');
            continue;
        end
    end
    
    stats = struct('meta',[],'index',[],'peaks',[]);
    try
        meta = cir_file.Strct_Metadata;
    catch me
        warning('problem with meta data  read')
        disp(me.message);
        continue;
    end
    
    % META DATA SECTION
    meta
    Ts = (1/meta.SampleRate_MHz_num)*1e-6;  % sample rate
    wl =  meta.CodewordLength_num;      % codeword length
    t = (0:Ts:Ts*(wl-1));% -Ts*pn_over;    % time array over a burst transmission
    t_view_max_ns = 1e3;            % maximum view in microseconds
    
    %NN = apf*rpa;
    NN = size(cir_file.IQdata,2);
    index = 1:NN;

    % setup output directories
    stats_dir = 'stats';
    fig_dir = 'figs';
    png_dir = 'png';
    mkdir(stats_dir);
    mkdir(fig_dir);
    mkdir(png_dir);

    % Compute some metrics for each CIR's
    ERP_dBW = 10*log10(meta.TransmitterPower_watts_num) ...
        + meta.TransmitterAntennaGain_dBi_num;
    ERP_W = 10^(ERP_dBW/10);
    ERP_Vrms = sqrt(ERP_W);
    Tx_rms_amp = sqrt(ERP_Vrms);
    for kk = 1:NN
        cir = cir_file.IQdata(:,kk);
        cir = [cir(1:end-20); zeros(20,1)]; 
        cir_mag = abs(cir);
        if ~isempty(cir_mag)
            if length(cir_mag) < wl
                peaks_t(kk) = NaN;
                peaks_k(kk) = NaN;
                K(kk) = NaN;   
                pdp_pwr(kk) = nf;  
                delay_spread(kk) = NaN;
                continue;
            end
            
            % compute the noise floor of the tail
            nf = mean(cir_mag(end-100:end-50));
            
            % compute the peak and the time of the peak
            [peaks(kk), peaks_k(kk)] = max(cir_mag(1:end-1024)); %#ok<*SAGROW>
            
            if peaks(kk) > nf*10
                
                % only select components 6 dB above the noise floor
                gtnf = cir_mag>nf*4;
                cir_gtnf = cir(gtnf);
                t_gtnf = t(gtnf);
                
                % compute the actual times of the peaks
                peaks_t(kk) = peaks_k(kk)*Ts;
                
                % compute the K factor assuming Rician channel
                r = cir_file.IQdata_Range_m(kk,3);
                K(kk) = compute_k_factor(t, cir, r, 6);
                
                % compute the duration of the CIR from the time of the
                % first component to time of the last component
                T = t_gtnf(end)-t_gtnf(1);
                Ngtnf = length(T);
                
                % compute the total power in the PDP
                pdp_pwr(kk) = sum(abs(Tx_rms_amp*cir_gtnf).^2)/Ngtnf; %#ok<*AGROW>
                
                % compute delay spread of the CIR
                delay_spread(kk) = compute_delay_spread(t_gtnf, cir_gtnf);
            else
                peaks_t(kk) = NaN;
                K(kk) = NaN;   
                pdp_pwr(kk) = nf;    
                delay_spread(kk) = NaN;
            end
        else
            peaks_t(kk) = NaN;
            K(kk) = NaN;
            pdp_pwr(kk) = nf;
            delay_spread(kk) = NaN;
        end
    end
    
    % Analyze the delay spread of the CIR's 
    if ~figvis, h = figure('Visible','on'); else h = figure(); end      
    [counts,centers] = hist(delay_spread*1e9,30);
    bar(centers, counts/sum(counts));
    xlabel('Delay Spread, Tau (nanosecs)')
    ylabel('Pr.(Tau)')
    grid on
    grid minor
    title('Histogram of Delay Spread')
    title({'Histogram of Delay Spread', strrep(mat_fname,'_','-')})
    drawnow
    savefig(h,[fig_dir '\' mat_fname(1:end-4) '__DelaySpread.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__DelaySpread.png'],'-dpng')    
    close(h)

    % Analyze the Rician K-factor 
    if ~figvis, h = figure('Visible','on'); else h = figure(); end      
    [counts,centers] = hist(K,30);
    bar(centers, counts/sum(counts));
    xlabel('K (dB)')
    ylabel('Pr.(K)')
    grid on
    grid minor
    title('Histogram of K-factor')
    title({'Histogram of K-factor', strrep(mat_fname,'_','-')})
    drawnow
    savefig(h,[fig_dir '\' mat_fname(1:end-4) '__Khist.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__Khist.png'],'-dpng')
    close(h)
    
    % Plot the CIR Magnitude
    if 0
    if ~figvis, h = figure('Visible','on'); else h = figure(); end 
    for kk = 1:NN
        cir = cir_file.IQdata(:,kk);
        if ~isempty(cir)
            plot(t*1e9, 10*log10(abs(cir)))
            title(sprintf('#%d CIR Mag', kk))
            xlim([0 t_view_max_ns]) 
            xlabel('time (ns)')
            ylabel('Gain (dB)')
            drawnow
        end
    end
    savefig(h, [fig_dir '\' mat_fname(1:end-4) '__cir_mag.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__cir_mag.png'],'-dpng')
    close(h)
    end

%     % View the time of peaks in time order
%     if ~figvis, h = figure('Visible','on'); else h = figure(); end
%     plot(peaks_t*1e9, 'd')
%     xlabel('record #')
%     ylabel('time (ns)')
%     title({'Time of Peak', strrep(mat_fname,'_','-')})
%     drawnow
%     savefig(h, [fig_dir '\' mat_fname(1:end-4) '_peaktime.fig']);
%     print(h,[png_dir '\' mat_fname(1:end-4) '__peak_time.png'],'-dpng')
%     close(h)

    % Analyzer receive power
    if ~figvis, h = figure('Visible','on'); else h = figure(); end
    plot(cir_file.IQdata_Range_m(:,3), 10*log10(pdp_pwr), 'o')
    xlabel('Distance (m)')
    ylabel('Received Power (dBW)')
    title({'Received Power', strrep(mat_fname,'_','-')})
    drawnow
    savefig(h, [fig_dir '\' mat_fname(1:end-4) '__power.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__power.png'],'-dpng')
    close(h)

    % save the metrics
    stats.meta = meta;
    stats.index = index;
    stats.pdp_pwr = pdp_pwr;
    stats.peaks = peaks; 
    stats.K = K; 
    stats.delay_spread = delay_spread;
    save([stats_dir '\' mat_fname(1:end-4) '__channel_stats.mat'], 'stats')

    % explicit clear of large memory
    cir_file = [];
    cir = [];
    cir_mag = [];
    K = [];
    gtnf = [];
    index = [];
    pdp_pwr = [];
    peaks = [];
    peaks_t = [];
    t = [];
    delay_spread = [];
    
    % save semaphore
    semv = 1;
    if ~exist('semv', 'dir')
        mkdir('semv');
    end
    save(sem_fname, 'semv')
    
    % gather memory stats
%     disp('memory usage')
%     memory

end

end % function
