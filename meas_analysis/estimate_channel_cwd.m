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
    % ANTENNA DATA
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
        % 10 dB above the noise floor.  If no samples exist, then we assume
        % that the cir is too noisy for channel estimation.
        %
        if any(peaks(kk) > nf*20)
            
            % only select components 6 dB above the noise floor
            gtnf = cir_mag>nf*4;
            cir_gtnf = cir(gtnf);
            t_gtnf = t(gtnf);

            % compute the actual times of the peaks
            peaks_t(kk) = peaks_k(kk)*Ts;

            % compute the K factor assuming Rician channel
            % we only consider components 10 dB above the noise floor and
            % where the peak occurs within 8 samples of beginning of the CIR 
            r = cir_file.IQdata_Range_m(kk,3);
            K(kk) = compute_k_factor(t, cir, r, 10, 8*Ts);

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
            % because of the wrapping of energy in the FFT-based
            % correlation we must remove the trailing edge.
            t_ds = t_gtnf(t_gtnf<(0.8*t_gtnf(end)));
            if length(t_ds) > 10  % we want at least 10 samples
                cir_ds = cir_gtnf(1:length(t_ds));
                [mean_delay_sec(kk), rms_delay_spread_sec(kk), cir_duration(kk)] = ...
                    compute_delay_spread(t_ds, cir_ds);
            end

        end
    end
    
    %
    % Extract range data for off-line analysis
    %
    r = cir_file.IQdata_Range_m(:,3);
    
    %
    % Analyze path loss versus distance
    %
    h = figure();
    k_gt3 = find(r>3);
    r_gt3 = r(k_gt3);
    r_p = r_gt3;  pl_p = path_gain_dB(k_gt3);
    r_p(isnan(pl_p)) = [];
    pl_p(isnan(pl_p)) = [];
    pl_p_fit = polyfit(r_p,pl_p,2);
    
    % one-segment linear fit
    p1=polyfit(r_p,pl_p,1); 
    r_px = linspace(min(r_p),max(r_p),100);
    semilogx(r_p, pl_p, 'o', r_px, polyval(p1, r_px));
    legend('Measured Data', ...
        sprintf('p=%0.2fx + %0.2f',p1));
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
    h = figure();      
    [counts,centers] = hist(1e9*rms_delay_spread_sec,50000);
    ds_probs = cumsum(counts/sum(counts));
    ds_centers = centers(ds_probs<0.99);
    ds_probs = ds_probs(ds_probs<0.99);
    plot(ds_centers, ds_probs);
    xlabel('Delay Spread, \tau_0 (nanosecs)')
    ylabel('Pr.(\tau < \tau_0)')
    setCommonGridProps()
    title({'CDF of Delay Spread', strrep(mat_fname,'_','-')})
    drawnow
    savefig(h,[fig_dir '\' mat_fname(1:end-4) '__ds.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__ds.png'],'-dpng')    
    close(h)
    
    %
    % Analyze the delay spread versus Distance
    %
    h = figure();      
    % remove nans from data
    k_gt3 = find(r>3);
    r_gt3 = r(k_gt3);
    r_p = r_gt3;  ds_p = rms_delay_spread_sec(k_gt3);
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
    % Analyze the Rician K Factor Estimates
    % 
    h = figure();      
    [counts,centers] = hist(K,30);
    plot(centers, cumsum(counts/sum(counts)));
    xlabel('K_0 (dB)')
    ylabel('Pr.(K < K_0)')
    setCommonGridProps()
    title('Histogram of K-factor')
    title({'CDF of K-factor', strrep(mat_fname,'_','-')})
    drawnow
    savefig(h,[fig_dir '\' mat_fname(1:end-4) '__Kcdf.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__Kcdf.png'],'-dpng')
    close(h)
    
    % 
    % Analyze the Rician K Factor Estimates Versus Distance
    % 
    h = figure();     
    r = cir_file.IQdata_Range_m(:,3);
    k_gt3 = find(r>3);
    r_gt3 = r(k_gt3);
    r_p = r_gt3;  K_p = K(k_gt3);
    r_p(isnan(K_p)) = [];
    K_p(isnan(K_p)) = [];
    plot(r_p, K_p, 'o');
    xlabel('distance (m))')
    ylabel('K (dB)')
    grid on
    grid minor
    title({'Rician K versus Distance', strrep(mat_fname,'_','-')})
    drawnow
    savefig(h,[fig_dir '\' mat_fname(1:end-4) '__KvRange.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__KvRange.png'],'-dpng')    
    close(h)
    
    % save the metrics
    stats = struct(...
        'meta',meta,...
        'path_gain_dB',path_gain_dB,...
        'path_gain_dB_poly',pl_p_fit,...
        'peaks',peaks,...
        'K',K,...
        'rms_delay_spread_sec',rms_delay_spread_sec, ...
        'mean_delay_sec',mean_delay_sec, ...
        'cir_duration_sec',cir_duration);
    
    save([stats_dir '\' mat_fname(1:end-4) '__channel_stats.mat'], 'stats')

    % explicit clear of large memory
    cir_file = []; %#ok<NASGU>
    
    % save semaphore
    semv = 1; %#ok<NASGU>
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



