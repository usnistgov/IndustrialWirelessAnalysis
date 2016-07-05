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
% query the list of measurement files
%
arr_dir = '.';          % sub directory of stored cir records
files = dir(pattern);

%
% Process each file in turn
%
for fk = 1:length(files)
    
    % use test data or the real thing
    try 
        
        mat_fname = files(fk).name; 
        cir_file_path = [arr_dir '\' mat_fname];       
        if TESTING
            disp('opening TEST_DATA');
            cir_file = TEST_DATA;
        else
            disp(['opening ' mat_fname]);
        
            % check for semphore file to indicate that it is already processed
            sem_fname = ['semv\' mat_fname(1:end-4) '__sem.mat'];
            if ~doall && exist(sem_fname, 'file')
                disp('semaphore file exists. overwriting...');
                cir_file = load(cir_file_path);
            elseif ~doall && exist(sem_fname, 'file')
                disp('semaphore file exists. skipping...');
                continue;
            else
                disp('semaphore file does not exist. loading...');
                cir_file = load(cir_file_path);
            end
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
    
    try
        meta = cir_file.Strct_Metadata;
    catch me
        warning('problem with meta data  read')
        disp(me.message);
        continue;
    end
    
    % META DATA SECTION
    disp(meta)
    Ts = (1/meta.SampleRate_MHz_num)*1e-6;  % sample rate
    wl =  meta.CodewordLength_num;      % codeword length
    t = (0:Ts:Ts*(wl-1));% -Ts*pn_over;    % time array over a burst transmission
    t_view_max_ns = 1e3;            % maximum view in microseconds
    
    %NN = apf*rpa;
    NN = size(cir_file.IQdata,2);
    index = 1:NN;

    % setup some memory for the metrics
    peaks = nan(NN,1);
    peaks_t = nan(NN,1);
    peaks_k = nan(NN,1);
    K = nan(NN,1);   
    path_gain_dB = nan(NN,1);  
    rms_delay_spread_sec = nan(NN,1);     
    mean_delay_sec = nan(NN,1); 
    cir_duration = nan(NN,1);    

    % setup output directories
    stats_dir = 'stats';
    fig_dir = 'figs';
    png_dir = 'png';
    mkdir(stats_dir);
    mkdir(fig_dir); 
    mkdir(png_dir); 

    %
    % META DATA
    % 
    TransmitterAntennaGain_dBi = meta.TransmitterAntennaGain_dBi_num;
    ReceiverAntennaGain_dBi = meta.ReceiverAntennaGain_dBi_num;
    
    % compute the ERP and RMS amplitude
    for kk = 1:NN
        
        % extract the CIR for this record from the data file
        cir = cir_file.IQdata(:,kk);
        
        % compute the magnitude of the CIR samples
        cir_mag = abs(cir);
        
        % ignore record if it is empty or the length is less than the
        % expected codeword length.  This indicates that something went
        % wrong with the instrumentation.
        if isempty(cir_mag)
            continue;
        elseif length(cir_mag) < wl
            continue;
        end

        % compute the noise floor of the tail end of the record but not the
        % mathematical wrapping of the forward impulse components 
        nf = mean(cir_mag(end-100:end-50));

        % compute the peak and the time of the peak of the CIR
        % save the peak information for later analysis
        [peaks(kk), peaks_k(kk)] = max(cir_mag(1:end-1024)); %#ok<*SAGROW>

        %
        % We only consider components CIR's where the peak is greater than
        % 6 dB above the noise floor.  If no samples exist, then we assume
        % that the cir is too noisy for channel estimation.
        %
        if any(peaks(kk) > nf*4)
            
            % only select components 6 dB above the noise floor
            gtnf = cir_mag>nf*4;
            cir_gtnf = cir(gtnf);
            t_gtnf = t(gtnf);

            % compute the actual times of the peaks
            peaks_t(kk) = peaks_k(kk)*Ts;

            % compute the K factor assuming Rician channel
            r = cir_file.IQdata_Range_m(kk,3);
            K(kk) = compute_k_factor(t, cir, r, 6);

            % Compute the path loss in the cir
            % note that the CIR contains antenna gains.  We must remove the
            % bulk antenna gains using the assumption that the gain is
            % applied equally to all multi-path components.  We know that
            % this is not the true case, but without ray-tracing is the
            % only option.
            Ngtnf = length(t_gtnf);
            path_gain_dB(kk) = 10*log10(sum(cir_gtnf.*conj(cir_gtnf))/Ngtnf) ...
                - TransmitterAntennaGain_dBi - ReceiverAntennaGain_dBi;  

            % compute delay spread parameters of the CIR 
            [mean_delay_sec(kk), rms_delay_spread_sec(kk), cir_duration(kk)] = ...
                compute_delay_spread(t_gtnf, cir_gtnf);

        end
    end
    
    %
    % Extract range data for off-line analysis
    %
    r = cir_file.IQdata_Range_m(:,3);
    
    %
    % Analyze path loss versus distance
    %
    if ~figvis, h = figure('Visible','on'); else h = figure(); end
    r_p = r;  pl_p = path_gain_dB;
    r_p(isnan(pl_p)) = [];
    pl_p(isnan(pl_p)) = [];
    semilogx(r_p, pl_p, 'o')
    setCommonGridProps()    
    xlabel('Distance (m)')
    ylabel('Path Loss (dB)')
    title({'Path Loss', strrep(mat_fname,'_','-')})
    drawnow
    savefig(h, [fig_dir '\' mat_fname(1:end-4) '__pathloss.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__pathloss.png'],'-dpng')
    close(h)    
    
    %
    % Analyze the delay spread of the CIR's 
    %
    if ~figvis, h = figure('Visible','on'); else h = figure(); end      
    [counts,centers] = hist(1e9*rms_delay_spread_sec,50000);
    ds_probs = cumsum(counts/sum(counts));
    ds_centers = centers(ds_probs<0.99);
    ds_probs = ds_probs(ds_probs<0.99);
    plot(ds_centers, ds_probs);
    xlabel('Delay Spread, Tau (nanosecs)')
    ylabel('Pr.(ds < \tau)')
    setCommonGridProps()
    title({'Cum Prob. of Delay Spread', strrep(mat_fname,'_','-')})
    drawnow
    savefig(h,[fig_dir '\' mat_fname(1:end-4) '__ds.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__ds.png'],'-dpng')    
    close(h)
    
    %
    % Analyze the delay spread versus Distance
    %
    if ~figvis, h = figure('Visible','on'); else h = figure(); end      
    % remove nans from data
    r_p = r;  ds_p = rms_delay_spread_sec;
    r_p(isnan(ds_p)) = [];
    ds_p(isnan(ds_p)) = [];
    plot(r_p,ds_p,'o');
    xlabel('distance (m))')
    ylabel('\tau (ns)')
    setCommonGridProps()
    title({'Delay Spread versus Distance', strrep(mat_fname,'_','-')})
    drawnow
    savefig(h,[fig_dir '\' mat_fname(1:end-4) '__ds2dist.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__ds2dist.png'],'-dpng')    
    close(h)
    
    %
    % Analyze the duration of the CIR's 
    %
    if ~figvis, h = figure('Visible','on'); else h = figure(); end      
    [counts,centers] = hist(1e9*cir_duration,1000);
    dur_probs = cumsum(counts/sum(counts));
    dur_centers = centers(dur_probs<0.99);
    dur_probs = dur_probs(dur_probs<0.99);
    plot(dur_centers, dur_probs);
    xlabel('cir duration, Tau (nanosecs)')
    ylabel('Pr.(dur < Tau)')
    setCommonGridProps()
    title({'Cum Prob. of CIR Duration', strrep(mat_fname,'_','-')})
    drawnow
    savefig(h,[fig_dir '\' mat_fname(1:end-4) '__dur.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__dur.png'],'-dpng')    
    close(h)    

    % 
    % Analyze the Rician K Factor Estimates
    % 
    if ~figvis, h = figure('Visible','on'); else h = figure(); end      
    [counts,centers] = hist(K,30);
    bar(centers, counts/sum(counts));
    xlabel('K (dB)')
    ylabel('Pr.(K)')
    setCommonGridProps()
    title('Histogram of K-factor')
    title({'Histogram of K-factor', strrep(mat_fname,'_','-')})
    drawnow
    savefig(h,[fig_dir '\' mat_fname(1:end-4) '__Khist.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__Khist.png'],'-dpng')
    close(h)
    
    % 
    % Analyze the Rician K Factor Estimates Versus Distance
    % 
    if ~figvis, h = figure('Visible','on'); else h = figure(); end     
    r = cir_file.IQdata_Range_m(:,3);
    % remove nans from data
    r_p = r;  K_p = K;
    r_p(isnan(K_p)) = [];
    K_p(isnan(K_p)) = [];
    plot(r_p, K_p, 'o');
    xlabel('distance (m))')
    ylabel('K (dB)')
    grid on
    grid minor
    title({'K versus Distance', strrep(mat_fname,'_','-')})
    drawnow
    savefig(h,[fig_dir '\' mat_fname(1:end-4) '__KvRange.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__KvRange.png'],'-dpng')    
    close(h)
    
    % 
    % Plot the CIR Magnitude
    % 
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
    
    % save the metrics
    stats = struct(...
        'meta',meta,...
        'path_gain_dB',path_gain_dB,...
        'peaks',peaks,...
        'K',K,...
        'rms_delay_spread_sec',rms_delay_spread_sec, ...
        'mean_delay_sec',mean_delay_sec, ...
        'cir_duration_sec',cir_duration);
    
    save([stats_dir '\' mat_fname(1:end-4) '__channel_stats.mat'], 'stats')

    % explicit clear of large memory
    cir_file = [];
    
    % save semaphore
    semv = 1;
    if ~exist('semv', 'dir')
        mkdir('semv');
    end
    save(sem_fname, 'semv')
    
    % gather memory stats
    % disp('memory usage')
    % memory

end

end % function

function setCommonGridProps
    grid on
    set(gca,'GridAlpha',0.5)
    set(gca,'MinorGridAlpha',0.5)
end



