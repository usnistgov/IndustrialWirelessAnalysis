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
OPT_AVGCIR_NTAP = 4;
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

if OPTS(OPT_AVGCIR_NTAP) == 1
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
        else
            if testForStatsFile(mat_fname(1:end-4));
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
    
    % META DATA SECTION
    disp(meta)
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
    path_gain_dB_poly = 0;
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
    mkdir(stats_dir);
    mkdir(fig_dir); 
    mkdir(png_dir); 

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
            if pk_pwr > -100;
                [mean_delay_sec(kk), rms_delay_spread_sec(kk), cir_duration(kk)] = ...
                    compute_delay_spread(Ts, cir);
            end
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
        if OPTS(OPT_KFACTOR) || OPTS(OPT_AVGCIR_NTAP)
            [K(kk), LOS(kk), k_pks] = compute_k_factor(cir, ns);
        end
        
        % Aggregate the sums for later computation of avg CIR
        % LOS and NLOS are considered as separate classes of CIR's
        if OPTS(OPT_AVGCIR_NTAP)
            if pk_pwr > -100
                cir0 = cir/max(abs(cir));
                if LOS(kk) == 1
                    num_los = num_los + 1;
                    inds = k_pks-k_pks(1)+1;
                    cir_sum_los(inds) = cir_sum_los(inds) + cir0(k_pks);
                elseif LOS(kk) == -1
                    num_nlos = num_nlos + 1;
                    inds = k_pks-k_pks(1)+1;
                    cir_sum_nlos(inds) = cir_sum_nlos(inds) + cir0(k_pks);
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
        
%         rA = cir_file.IQdata_Range_m(1:NN,3);
%         r = rA;
%         path_gain_dB = path_gain_dB(r~=r(1));
%         r = r(r~=r(1));
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
        
        infpt = 30;
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

        % upper segment
        Gfit2 = G(R>infpt);
        Rfit2 = R(R>infpt);
        p2 = polyfit(log10(Rfit2), Gfit2, 1);   
        path_gain_dB_poly = p2;

        % find draw points
        if ~isempty(Rfit1)
            x_intersect = fzero(@(x) polyval(p1-p2,x),log10(infpt));
            %x_intersect = fzero(@(x) polyval(p1-p2,x),[log10(min(Rfit1)), log10(max(Rfit2))]);
            d_val1 = log10(logspace(log10(min(Rfit1)), x_intersect, 20));
            g_val1 = polyval(p1, d_val1);
        else
            x_intersect = log10(min(R));
        end
        d_val2 = log10(logspace(x_intersect, log10(max(Rfit2)), 20));
        g_val2 = polyval(p2, d_val2);
        
        % frii as reference
        d_frii = logspace(log10(min(R)), log10(max(R)), 20);
        ff = meta.Frequency_GHz_num;
        c = physconst('LightSpeed');
        frii_fspl_dB = 10*log10(d_frii.^2) + 20*log10(ff) + 20*log10(1e9) + 20*log10(4*pi/c);         
        
        % plot the gains
        h = figure();
        semilogx(R,G, 'color', [0,0,0]+0.7, 'marker', '.', 'linestyle' , 'none')
        hold on
        if ~isempty(Rfit1)
            semilogx(10.^d_val1,g_val1,'b+-')
        end
        semilogx(10.^d_val2,g_val2,'b+-')
        semilogx(d_frii, -frii_fspl_dB, 'b-')
        text(d_frii(floor(length(d_frii)/4)),-frii_fspl_dB(floor(length(d_frii)/4))+5,'FSPL','Color','blue');         
        hold off
        
        setCommonAxisProps()    
        xlabel('distance (m)')
        ylabel('Path Gain (dB)')
        drawnow
        savefig(h, [fig_dir '\' mat_fname(1:end-4) '__pl.fig']);
        setFigureForPrinting();
        print(h,[png_dir '\' mat_fname(1:end-4) '__pl.png'],'-dpng','-r300')
        close(h) 

%         path_gain_calc = path_gain_dB;
%         r_p = r;  pl_p = path_gain_calc;
%         r_p = r_p(~isnan(pl_p));
%         pl_p = pl_p(~isnan(pl_p));
%         r_min_fit = 3;
%         r_max_fit = 0.9*max(r_p);    
%         k_lfit = find(r_p>r_min_fit & r_p<r_max_fit);
%         r_p_fit = r_p(k_lfit);
%         pl_p_fit = pl_p(k_lfit);
%         path_gain_dB_poly = polyfit(log10(r_p_fit),pl_p_fit,1);
%         this_path_gain_poly = path_gain_dB_poly;
% 
%         h = figure();
%         stdPathGain = std(path_gain_calc(~isnan(path_gain_calc)));
%         r_p_plot = log10(logspace(log10(min(r_p_fit)), log10(max(r_p_fit)), 10));
%         pl_poly_vals = polyval(this_path_gain_poly, r_p_plot);
%         semilogx(r_p, pl_p, 'color', [0,0,0]+0.7, 'marker', '.', 'linestyle' , 'none'); 
%         hold on
%         semilogx(10.^r_p_plot, pl_poly_vals, 'k-', ...
%            10.^r_p_plot, repmat(pl_poly_vals(:),1,2)+stdPathGain*[ones(10,1) -ones(10,1)], ...
%             'k--', 'LineWidth', 1.0);
%         
%         % frii as reference
%         if 1
%         d = r_p_plot;
%         ff = meta.Frequency_GHz_num;
%         c = physconst('LightSpeed');
%         frii_fspl_dB = 10*log10(d.^2) + 20*log10(ff) + 20*log10(1e9) + 20*log10(4*pi/c);
%         semilogx(d, -frii_fspl_dB, 'b-')
%         text(d(floor(length(d)/4)),-frii_fspl_dB(floor(length(d)/4))+5,'FSPL','Color','blue');
%         end       
%         
%         hold off
%         if 0
%         legend({'measured', ...
%             sprintf('%0.2fx + %0.1f',path_gain_dB_poly), ...
%             '+/- \sigma'}, 'Location', 'best');
%         end
%         setCommonAxisProps()    
%         xlabel('distance (m)')
%         ylabel('Path Gain (dB)')
%         drawnow
%         savefig(h, [fig_dir '\' mat_fname(1:end-4) '__pl.fig']);
%         setFigureForPrinting();
%         print(h,[png_dir '\' mat_fname(1:end-4) '__pl.png'],'-dpng','-r300')
%         close(h)  
        
    end % if OPTS(OPT_PATH_GAIN)
    

    if OPTS(OPT_DELAY_SPREAD) && sum(~isnan(mean_delay_sec)) > 100
        %
        % Analyze the average delay of the CIR's 
        %
        h = figure();      
        ds_th = nanmean(mean_delay_sec) + 3*nanstd(mean_delay_sec);
        [du_counts,du_centers] = hist(1e9*mean_delay_sec(mean_delay_sec<ds_th),50);
        du_probs = cumsum(du_counts/sum(du_counts));
        du_centers = du_centers(du_probs < 0.99);
        du_probs = du_probs(du_probs < 0.99);        
        yyaxis left, area(du_centers, [0 diff(du_probs)],'FaceAlpha',0.5)
        str = 'average delay, $$\tau_D$$ (ns)';xlabel(str,'Interpreter','Latex')
        yyaxis right, plot(du_centers,du_probs)
        str = 'Pr. $$ \hat{\tau_D} < {\tau_D} $$'; ylabel(str,'Interpreter','Latex'); 
        str = 'mean delay, $${\tau_D}$$ (ns)';xlabel(str,'Interpreter','Latex')
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
        ds_th = nanmean(rms_delay_spread_sec) + 3*nanstd(rms_delay_spread_sec);
        [ds_counts,ds_centers] = hist(1e9*rms_delay_spread_sec(rms_delay_spread_sec<ds_th),50);
        ds_probs = cumsum(ds_counts/sum(ds_counts));
        ds_centers = ds_centers(ds_probs < 0.99);
        ds_probs = ds_probs(ds_probs < 0.99);
        yyaxis left, area(ds_centers, [0 diff(ds_probs)],'FaceAlpha',0.5)
        str = 'Pr. $$\hat{S}$$'; ylabel(str,'Interpreter','Latex');
        yyaxis right, plot(ds_centers,ds_probs)
        str = 'Pr. $$\hat{S} < S$$'; ylabel(str,'Interpreter','Latex');
        str = 'rms delay spread, $$S$$ (ns)';xlabel(str,'Interpreter','Latex')
        setCommonAxisProps();
        drawnow
        savefig(h,[fig_dir '\' mat_fname(1:end-4) '__ds.fig']);
        h.PaperPositionMode = 'auto';
        setFigureForPrinting();
        print(h,[png_dir '\' mat_fname(1:end-4) '__ds.png'],'-dpng','-r300')    
        close(h)

        %
        % Analyze the delay spread versus distance
        %
        if 0
        h = figure();      
        % remove nans from data
        k_lfit = find(r>10);
        r_lfit = r(k_lfit);
        r_p = r_lfit;  
        ds_p = rms_delay_spread_sec(k_lfit);
        r_p(isnan(ds_p)) = [];
        ds_p(isnan(ds_p)) = [];
        plot(r_p,1e9*ds_p, 'color', [0,0,0]+0.7, ...
            'marker', '.', 'linestyle' , 'none');
        hold on;
        stdS = 1e9*std(rms_delay_spread_sec(~isnan(rms_delay_spread_sec)));
        ds_poly = polyfit(r_p,1e9*ds_p,1);
        r_p_plot = linspace(min(r_p), max(r_p), 10);
        ds_poly_vals = polyval(ds_poly, r_p_plot);
        plot(r_p_plot, ds_poly_vals, 'k-', ...
           r_p_plot, repmat(ds_poly_vals(:),1,2)+stdS*[ones(10,1) -ones(10,1)], ...
            'k--', 'LineWidth', 1.0);
        hold off;
        setCommonAxisProps()
        legend({'measured', ...
            sprintf('%0.2fx + %0.1f',ds_poly), ...
            '+/- \sigma'}, 'Location', 'best');        
        xlabel('distance, d (m)')
        ylabel('S (ns)')
        grid off
        %title({'RMS Delay Spread versus Distance', strrep(mat_fname,'_','-')})
        drawnow
        savefig(h,[fig_dir '\' mat_fname(1:end-4) '__ds2dist.fig']);
        setFigureForPrinting();
        print(h,[png_dir '\' mat_fname(1:end-4) '__ds2dist.png'],'-dpng','-r300')    
        close(h)  
        end
        
    end %if OPTS(OPT_DELAY_SPREAD)

    if OPTS(OPT_KFACTOR)
        % 
        % Analyze the Rician K Factor Estimates CDF
        % 
        h = figure();      
        Kx = K(~isinf(K));
        [counts,centers] = hist(Kx,300);
        plot(centers, cumsum(counts/sum(counts)), 'k');
        ylim([0 1]);
        str = '$$K$$ (dB)'; xlabel(str,'Interpreter','Latex');
        str = '$$ \hat{K} < K $$'; ylabel(str,'Interpreter','Latex');
        setCommonAxisProps()
        %title({'CDF of Rician K-factor Estimate', strrep(mat_fname,'_','-')})
        drawnow
        savefig(h,[fig_dir '\' mat_fname(1:end-4) '__Kcdf.fig']);
        setFigureForPrinting();
        print(h,[png_dir '\' mat_fname(1:end-4) '__Kcdf.png'],'-dpng','-r300')
        close(h)

        % 
        % Analyze the Rician K Factor Estimates Versus Distance
        % 
        if 0
        h = figure();     
        r = rA;
        k_lfit = find(r>3);
        r_lfit = r(k_lfit);
        r_p = r_lfit;  K_p = K(k_lfit);
        not_these = isnan(K_p) | isinf(K_p);
        r_p(not_these) = [];
        K_p(not_these) = [];
        plot(r_p, K_p, 'color', [0,0,0]+0.7, 'marker', '.', 'linestyle' , 'none');
        hold on 
        KdB_poly = polyfit(r_p,K_p,1);
        r_p_plot = linspace(min(r_p), max(r_p), 10);
        KdB_poly_vals = polyval(KdB_poly, r_p_plot);
        stdK = std(K_p);
        plot(r_p_plot, KdB_poly_vals, 'k-', ...
           r_p_plot, repmat(KdB_poly_vals(:),1,2)+2*stdK*[ones(10,1) -ones(10,1)], ...
            'k--', 'LineWidth', 1.0);
        hold off
        setCommonAxisProps()
        legend({'measured', ...
            sprintf('%0.2fx + %0.1f',KdB_poly), ...
        '+/- 2\sigma'},'FontSize',9);
        xlabel('distance (m)')
        ylabel('K (dB)')
        %title({'Rician K versus Distance', strrep(mat_fname,'_','-')})
        drawnow
        savefig(h,[fig_dir '\' mat_fname(1:end-4) '__KvRange.fig']);
        setFigureForPrinting();
        print(h,[png_dir '\' mat_fname(1:end-4) '__KvRange.png'],'-dpng','-r300')    
        close(h)   
        end
    
    end % OPTS(OPT_KFACTOR)
    
    % approximate an N-tap CIR from the measured CIR's
    if OPTS(OPT_AVGCIR_NTAP)
        
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
            %subplot(4,1,1:2)
            if OPTS(OPT_NTAP_APPROX)
                plot(1E9*t_ciravg, abs(cir_avg)); 
            else
                stem(1E9*t_ciravg, abs(cir_avg)); 
            end
            str = '$$\mid{h(t)}\mid$$';ylabel(str, 'Interpreter', 'Latex')
            xlabel('time (ns)')
            %set(gca,'XTickLabel','')
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
            
            if 0
            subplot(4,1,3:4)
            plot(1E9*t_ciravg, angle(cir_avg));
            str = '$$\angle{h(t)}$$';ylabel(str, 'Interpreter', 'Latex')
            xlabel('time (ns)')
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
            close(h)   
        
        end

        
    end % OPTS(OPT_AVGCIR_NTAP)
    
    if OPTS(OPT_WRITE_STATS)
        % save the metrics
        stats = struct(...
            'meta',meta,...
            'path_gain_range_m', path_gain_range_m, ...
            'path_gain_dB',path_gain_dB,...
            'path_gain_dB_poly',path_gain_dB_poly,...
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

function setCommonAxisProps()

    alw = 0.75;    % AxesLineWidth
    fsz = 10;      % Fontsize
    lw = 0.75;%1.25;      % LineWidth
    msz = 6.0;       % MarkerSize
    
%    grid on
    set(gca,'XGrid','on')
    set(gca,'XMinorGrid','off')
    set(gca,'YGrid','on')
    set(gca,'YMinorGrid','off')
    set(gca,'GridAlpha',0.25)
    set(gca,'MinorGridAlpha',0.4)
    set(gca,'Fontsize',fsz)
    set(gca,'LineWidth',alw);
    set(gca,'FontName','TimesRoman')
    
    % set the line properties
    hline = get(gca,'Children');
    for h = hline(:)'
        if strcmp(h.Type,'line')
            h.LineWidth = lw;
            h.MarkerSize = msz;
        end
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


