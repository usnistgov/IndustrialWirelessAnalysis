function estimate_channel_cwd(pattern, doall, figvis, Ts, TEST_DATA)
% Analyze complex impulse responses from measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

TESTING = false;
if nargin == 5
    TESTING = true;
end

if nargin < 4
    Ts = 12.5e-9;
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
    
    stats = struct('meta',[],'index',[],'peaks',[],'vars',[]);
    try
        meta = cir_file.Strct_Metadata;
    catch me
        warning('problem with meta data  read')
        disp(me.message);
        continue;
    end
    
    % META DATA SECTION
    apf = meta.NumberAcqusitions_num;       % acquisitions per file
    pn_over = meta.PNOversample_num;        % samples per chip
    rpa = meta.NumberRecordperAcqusition_num;      % records per acquisition
    ns =  meta.CodewordLength_num;      % number of samples
%     Ts = meta{14,2};                % sample period
    wl = ns;        % codeword length
    t = (0:Ts:Ts*(wl-1))-Ts*pn_over;    % time array over a burst transmission
    t_view_max_us = 1.6;            % maximum view in microseconds
    
    NN = apf*rpa;
    index = 1:NN;

    % setup output directories
    stats_dir = 'stats';
    fig_dir = 'figs';
    png_dir = 'png';
    mkdir(stats_dir);
    mkdir(fig_dir);
    mkdir(png_dir);

    % Compute some metrics for each CIR's
    Tx_pwr = NaN; %#ok<NASGU> % Watts
    if ~isempty(strfind(mat_fname, '2G'))
        Tx_pwr = 1.5;
    else % 5 GHz
        Tx_pwr = 1.25;
    end
    Tx_rms_amp = sqrt(Tx_pwr);
    for kk = 1:NN
        cir = cir_file.IQdata(:,kk);
        cir = circshift(cir, pn_over);
        cir = [cir(1:end-20); zeros(20,1)]; 
        cir_mag = abs(cir);
        if ~isempty(cir_mag)
            if length(cir_mag) < ns
                peaks_t(kk) = NaN;
                vars(kk) = NaN;
                K(kk) = NaN;   
                pdp_pwr(kk) = nf;            
                continue;
            end
            nf = mean(cir_mag(end-100:end-50));
            [peaks(kk), peaks_t(kk)] = max(cir_mag(1:end-1024)); %#ok<*SAGROW>
            if peaks(kk) > nf*10
                gtnf = cir_mag>nf*10;
                cir_gtnf = cir(gtnf);
                t_gtnf = t(gtnf);
                peaks_t(kk) = peaks_t(kk)*Ts;
                var_kk = var(cir_mag(cir_mag>(nf*4)));
                K(kk) = 10*log10(peaks(kk)^2/(2*var_kk));   
                T = t_gtnf(end)-t_gtnf(1);
                Ngtnf = length(T);
    %             pdp_pwr(kk) = sum(abs(Tx_rms_amp*cir_gtnf).^2)/T; %#ok<*AGROW>
                pdp_pwr(kk) = sum(abs(Tx_rms_amp*cir_gtnf).^2)/Ngtnf; %#ok<*AGROW>
            else
                peaks_t(kk) = NaN;
                vars(kk) = NaN;
                K(kk) = NaN;   
                pdp_pwr(kk) = nf;            
            end
        else
            peaks_t(kk) = NaN;
            vars(kk) = NaN;
            K(kk) = NaN;
            pdp_pwr(kk) = nf;
        end
    end

    % Analyze the Rician K-factor 
    if ~figvis, h = figure('Visible','off'); else h = figure(); end      
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
    
    

    % Plot the CIR Magnitude
    if 0
    if ~figvis, h = figure('Visible','off'); else h = figure(); end 
    for kk = index
        cir = cir_file.the_run.cir(:,kk);
        cir = circshift(cir, pn_over);
        cir = [cir(1:end-20); zeros(20,1)]; 
        if ~isempty(cir)
            plot(t*1e6, 10*log10(abs(cir)))
            title(sprintf('#%d CIR Mag', kk))
            xlim([0 t_view_max_us]) 
            xlabel('time (us)')
            ylabel('Gain (dB)')
            drawnow
        end
    end
    savefig(h, [fig_dir '\' mat_fname(1:end-4) '__cir_mag.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__cir_mag.png'],'-dpng')
    end

    % View the time of peaks in time order
    if ~figvis, h = figure('Visible','off'); else h = figure(); end
    plot(peaks_t*1e9, 'd')
    xlabel('record #')
    ylabel('time (ns)')
    title({'Time of Peak', strrep(mat_fname,'_','-')})
    drawnow
    savefig(h, [fig_dir '\' mat_fname(1:end-4) '_peaktime.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__peak_time.png'],'-dpng')

    % View the time of peaks in ascending distance
    % if ~figvis, h = figure('Visible','off'); else h = figure(); end
    % plot

    % view the K-factor as a function of ascending distance to transmitter
    % if ~figvis, h = figure('Visible','off'); else h = figure(); end
    % plot

    % Analyzer receive power
    if ~figvis, h = figure('Visible','off'); else h = figure(); end
    plot(index, 10*log10(pdp_pwr), 'o')
    xlabel('index')
    ylabel('Received Power (dBm)')
    title({'Received Power', strrep(mat_fname,'_','-')})
    drawnow
    savefig(h, [fig_dir '\' mat_fname(1:end-4) '__power.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__power.png'],'-dpng')

    % save the metrics
    stats.meta = meta;
    stats.index = index;
    stats.pdp_pwr = pdp_pwr;
    stats.peaks = peaks; 
    stats.K = K; %#ok<STRNU>
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
