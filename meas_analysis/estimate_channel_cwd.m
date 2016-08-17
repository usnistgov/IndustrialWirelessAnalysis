function estimate_channel_cwd(pattern, doall, TEST_DATA)
% Analyze complex impulse responses from measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

TESTING = false;
if nargin == 3
    TESTING = true;
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
    wl =  meta.CodewordLength_num;          % codeword length
    ns = meta.PNOversample_num;             % oversample rate
    t = (0:Ts:Ts*(wl-1));% -Ts*pn_over;     % time array over a burst transmission
    
    %NN = apf*rpa;
    NN = size(cir_file.IQdata,2);

    % setup some memory for the metrics
    peaks = nan(NN,1);
    peaks_t = nan(NN,1);
    peaks_k = nan(NN,1);
    K = nan(NN,1);  
    LOS = nan(NN,1);  
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
        
        % range of the measurement
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
        k = select_cir_samples(r, t, cir);
        t_k = t(k);
        cir_k = cir(k);

        % compute the noise floor of the tail end of the record but not the
        % mathematical wrapping of the forward impulse components 
        nf = mean(cir_mag2(end-100:end-50));

        % compute the peak and the time of the peak of the CIR
        % save the peak information for later analysis
        [peaks(kk), peaks_k(kk)] = max(cir_mag2(1:end-1024)); %#ok<*SAGROW>

        % compute the actual times of the peaks
        peaks_t(kk) = peaks_k(kk)*Ts;

        % Compute the path loss in the cir
        % note that the CIR contains antenna gains.  We must remove the
        % bulk antenna gains using the assumption that the gain is
        % applied equally to all multi-path components.  We know that
        % this is not the true case, but without ray-tracing it is the
        % only option.
        path_gain_dB(kk) = compute_path_gain(cir_k, ...
            TransmitterAntennaGain_dBi, ...
            ReceiverAntennaGain_dBi);           

        % compute the K factor assuming Rician channel
        % we only consider components 10 dB above the noise floor and
        % where the peak occurs within 8 samples of beginning of the CIR 
        [K(kk), LOS(kk)] = compute_k_factor(t_k, cir_k, ns);

        % compute delay spread parameters of the CIR 
        % because of the wrapping of energy in the FFT-based
        % correlation we must remove the trailing edge.
        [mean_delay_sec(kk), rms_delay_spread_sec(kk), cir_duration(kk)] = ...
            compute_delay_spread(t_k, cir_k, nf);
        if rms_delay_spread_sec(kk) >= 0.5e-5
            %disp(max(t_k)/mean(t_k));
            h=figure();
            plot(t,10*log10(abs(cir).^2/max(abs(cir).^2)),...
                t_k,10*log10(abs(cir_k).^2/max(abs(cir_k).^2)),'ro')
            refline(0, 10*log10(nf/max(abs(cir).^2)))
            close(h)
        end

    end
    
    %
    % Extract range data for off-line analysis
    %
    r = cir_file.IQdata_Range_m(:,3);
    path_gain_dB = path_gain_dB(r~=r(1));
    r = r(r~=r(1));
    
    %
    % Analyze path loss versus distance
    %
    r_p = r;  pl_p = path_gain_dB;
    r_p = r_p(~isnan(pl_p));
    pl_p = pl_p(~isnan(pl_p));
    r_min_fit = 1.05*min(r_p);
    r_max_fit = 0.95*max(r_p);    
    k_fit = find(r_p>r_min_fit & r_p<r_max_fit);
    r_p_fit = r_p(k_fit);
    pl_p_fit = pl_p(k_fit);
    p1 = polyfit(r_p_fit,pl_p_fit,1);
    
    h = figure();
    r_px = linspace(min(r_p),max(r_p),100);
    plot(r_p, pl_p, 'bo', r_px, polyval(p1, r_px), 'r-');
    legend('Measured Data', ...
        sprintf('p=%0.2fx + %0.2f',p1));
    setCommonAxisProps()    
    xlabel('distance [m]')
    ylabel('Path Gain (dB)')
    title({'Path Loss', strrep(mat_fname,'_','-')})
    drawnow
    savefig(h, [fig_dir '\' mat_fname(1:end-4) '__pl.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__pl.png'],'-dpng')
    close(h)    
    
    % path loss error bars
    h = figure();
    [r_bins, xx_u, xx_s] = makeErrorBars(gca(), r_p, pl_p, 20, 'log');
    set(gca,'xscale','log');
    setCommonAxisProps();
    title({'Path Loss (Error Bars)', strrep(mat_fname,'_','-')})
    xlabel('distance [m]')
    ylabel('Path Gain [dB]')
    drawnow
    savefig(h,[fig_dir '\' mat_fname(1:end-4) '__pl_eb.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__pl_eb.png'],'-dpng')    
    close(h)   
    
    % path loss log normal distribution, variance versus distance
    h = figure();
    stem(r_bins, xx_s);
    set(gca,'xscale','log');
    setCommonAxisProps();
    title({'Path Gain Variation versus Distance', strrep(mat_fname,'_','-')})
    xlabel('distance [m]')
    ylabel('\sigma [dB]')
    drawnow
    savefig(h,[fig_dir '\' mat_fname(1:end-4) '__plstd.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__plstd.png'],'-dpng')    
    close(h)       

    
    %
    % Analyze the average delay of the CIR's 
    %
    h = figure();      
    [du_counts,du_centers] = hist(1e9*mean_delay_sec,50000);
    du_probs = cumsum(du_counts/sum(du_counts));
    plot(du_centers, du_probs);
    ylim([0 1]);
    xlim([0 max(du_centers(du_probs<0.995))]);
    xlabel('average delay spread, t_D [ns]')
    ylabel('Pr.(T_D < t_D)')
    setCommonAxisProps();
    title({'CDF of Average Delay', strrep(mat_fname,'_','-')})
    drawnow
    savefig(h,[fig_dir '\' mat_fname(1:end-4) '__du.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__du.png'],'-dpng')    
    close(h)
    
    %
    % Analyze the delay spread of the CIR's 
    %
    h = figure();      
    [ds_counts,ds_centers] = hist(1e9*rms_delay_spread_sec,5000);
    ds_probs = cumsum(ds_counts/sum(ds_counts));
    plot(ds_centers, ds_probs);
    ylim([0 1]);
    xlabel('rms delay spread, s [ns]')
    ylabel('Pr.(S < s)')
    setCommonAxisProps();
    title({'CDF of Delay Spread', strrep(mat_fname,'_','-')})
    drawnow
    savefig(h,[fig_dir '\' mat_fname(1:end-4) '__ds.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__ds.png'],'-dpng')    
    close(h)
    
    %
    % Analyze the delay spread versus distance
    %
    h = figure();      
    % remove nans from data
    k_gt3 = find(r>3);
    r_gt3 = r(k_gt3);
    r_p = r_gt3;  
    ds_p = rms_delay_spread_sec(k_gt3);
    r_p(isnan(ds_p)) = [];
    ds_p(isnan(ds_p)) = [];
    plot(r_p,1e9*ds_p,'o');
    xlabel('distance, d [m])')
    ylabel('S [ns]')
    setCommonAxisProps()
    title({'RMS Delay Spread versus Distance', strrep(mat_fname,'_','-')})
    drawnow
    savefig(h,[fig_dir '\' mat_fname(1:end-4) '__ds2dist.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__ds2dist.png'],'-dpng')    
    close(h)   
    
    % delay spread error bars
    h = figure();  
    makeErrorBars(gca(), r_p, ds_p, 30, 'linear', 1e9);
    xlabel('distance, d [m])')
    ylabel('S [ns]')    
    setCommonAxisProps()
    title({'RMS Delay Spread versus Distance (Error Bars)', strrep(mat_fname,'_','-')})
    drawnow
    savefig(h,[fig_dir '\' mat_fname(1:end-4) '__ds2dist_eb.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__ds2dist_eb.png'],'-dpng')    
    close(h)   

    % 
    % Analyze the Rician K Factor Estimates
    % 
    h = figure();      
    [counts,centers] = hist(K,30);
    plot(centers, cumsum(counts/sum(counts)));
    ylim([0 1]);
    xlabel('K_0 [dB]')
    ylabel('Pr.(K < K_0)')
    setCommonAxisProps()
    title({'CDF of Rician K-factor Estimate', strrep(mat_fname,'_','-')})
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
    xlabel('distance [m])')
    ylabel('K [dB]')
    grid on
    grid minor
    title({'Rician K versus Distance', strrep(mat_fname,'_','-')})
    drawnow
    savefig(h,[fig_dir '\' mat_fname(1:end-4) '__KvRange.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__KvRange.png'],'-dpng')    
    close(h)
    
    % Rician K error bars
    h = figure();  
    makeErrorBars(gca(), r_p, K_p, 25);
    setCommonAxisProps();
    xlabel('distance [m])')
    ylabel('K [dB]')    
    title({'Rician K (Error Bars)', strrep(mat_fname,'_','-')})
    drawnow
    savefig(h,[fig_dir '\' mat_fname(1:end-4) '__KvRange_eb.fig']);
    print(h,[png_dir '\' mat_fname(1:end-4) '__KvRange_eb.png'],'-dpng')    
    close(h)       
    
    % save the metrics
    stats = struct(...
        'meta',meta,...
        'path_gain_dB',path_gain_dB,...
        'path_gain_dB_poly',p1,...
        'peaks',peaks,...
        'K',K,...
        'los', LOS, ...
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

cmp_pl_poly( '*_stats.mat', '.\stats', '.\figs', '.\png' )

end % function

function setCommonAxisProps()
    grid on
    set(gca,'GridAlpha',0.5)
    set(gca,'MinorGridAlpha',0.5)
    set(gca,'Fontsize',12)
    set(gca,'FontName','TimesRoman')
end



function [r_bins, xx_u, xx_s] = makeErrorBars(ha, r, v, n, xscale, yscale)
    if nargin < 6
        yscale = 1;
    end
    if nargin < 5
        xscale = 'linear';
    end
    if strcmp(xscale, 'log')
        r_bins=logspace(log10(min(r)),log10(max(r)),n);
    else
        r_bins=linspace(min(r),max(r),n);
    end
    r_d = discretize(r,r_bins);
    Nbins = length(r_bins);
    xx_u = zeros(Nbins,1);
    xx_s = zeros(Nbins,1);
    for ii = 1:Nbins
        xx = v(r_d==ii);
        xx_u(ii) = mean(xx);
        xx_s(ii) = std(xx);
    end
    errorbar(ha, r_bins, yscale*xx_u, yscale*xx_s,'kd-');
end

