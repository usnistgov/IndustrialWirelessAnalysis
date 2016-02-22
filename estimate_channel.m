function estimate_channel(pattern)
% Analyze complex impulse responses from measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

more off; 
grid minor;
set(0,'DefaultFigureWindowStyle','docked')

files = dir(pattern);
for fk = 1:length(files)
    
mat_fname = files(fk).name;
disp(['opening ' mat_fname]);
cir_file = load(mat_fname);

%
peaks = [];
peaks_t = [];
K = [];

stats = struct('meta',[],'index',[],'peaks',[],'vars',[]);
try
    meta = cir_file.AllthePdps(1:19,7:8);
catch me
    disp me;
    continue;
end
apf = meta{11,2};
pn_over = cell2mat(meta(7,2));
rpa = cell2mat(meta(9,2));
ns =  cell2mat(meta(8,2));
Ts = meta{14,2};
wl = meta{8,2};
t = [0:Ts:Ts*(wl-1)]-Ts*pn_over;
t_view_max_us = 1.6;

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
    cir = cell2mat(cir_file.AllthePdps(kk,6));
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
            
%             pdp_pwr(kk) = sum(abs(Tx_rms_amp*cir_gtnf).^2)/T; %#ok<*AGROW>
            pdp_pwr(kk) = sum(abs(Tx_rms_amp*cir_gtnf).^2); %#ok<*AGROW>
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

% Analyzer the Rician K-factor
h = figure(1); 
[counts,centers] = hist(K,30);
bar(centers, counts/sum(counts));
xlabel('K (dB)')
ylabel('Pr.(K)')
grid on
grid minor
title('Histogram of K-factor')
drawnow
savefig(h,[fig_dir '\' 'hist_' mat_fname '.fig']);
print(h,[png_dir '\' mat_fname '__Khist.png'],'-dpng')

% Plot the CIR Magnitude
if 0
h = figure(2); 
for kk = index
    cir = cell2mat(cir_file.AllthePdps(kk,6));
    cir = circshift(cir, pn_over);
    cir = [cir(1:end-20); zeros(20,1)]; 
    if ~isempty(cir)
        plot(t*1e6, 10*log10(abs(cir)))
        str_rectime = cell2mat(cir_file.AllthePdps(kk,5));
        title(sprintf('#%d CIR Mag at %s', kk, str_rectime))
        xlim([0 t_view_max_us]) 
        xlabel('time (us)')
        ylabel('Gain (dB)')
        drawnow
    end
end
savefig(h, [fig_dir '\' 'cirmag_' mat_fname '.fig']);
print(h,[png_dir '\' mat_fname '__cir_mag.png'],'-dpng')
end

% View the time of peaks in time order
h = figure(3); clf
plot(peaks_t*1e9, 'd')
xlabel('record #')
ylabel('time (ns)')
title('Time of Peak')
drawnow
savefig(h, [fig_dir '\' 'peaktime_' mat_fname '.fig']);
print(h,[png_dir '\' mat_fname '__peak_time.png'],'-dpng')

% View the time of peaks in ascending distance
% figure(4)
% plot

% view the K-factor as a function of ascending distance to transmitter
% figure(5)
% plot

% Analyzer receive power
h = figure(6);
plot(index, 10*log10(pdp_pwr), 'o')
xlabel('index')
ylabel('Received Power (dBm)')
title(strrep(mat_fname,'_','-'))
drawnow
savefig(h, [fig_dir '\' 'power_' mat_fname '.fig']);
print(h,[png_dir '\' mat_fname '__power.png'],'-dpng')

% save the metrics
stats.meta = meta;
stats.index = index;
stats.pdp_pwr = pdp_pwr;
stats.peaks = peaks; 
stats.K = K; %#ok<STRNU>
save([stats_dir '\' mat_fname '__channel_stats.mat'], 'stats')

close all;
clear('cir*')

end

end % function
