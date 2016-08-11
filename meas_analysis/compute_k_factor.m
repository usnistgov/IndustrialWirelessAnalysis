function [ K ] = compute_k_factor( t, cir, r, ns )
%COMPUTE_K_FACTOR Computes the K-factor of the impulse response
%   Computes the K factor of the CIR using all components of magnitude
%   greater than the noise floor times L.  The range, r, in meters guides
%   the extraction of the LOS power component.
%
%   t:      time in seconds
%   cir:    the channel impulse response
%   r:      range in meters (not used)
%   ns:    oversample rate
% 
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

Ts = t(2)-t(1);
dT0 = ns*Ts;

% eliminate the anomalous tail components (last 8 samples)
t = t(1:end-8);
cir = cir(1:end-8);

% select peak values above the noise floor and within 12 dB (15.8489) of
% the peak magnitude per ITU P1407-5 S2.2.7
% the PDP is normalized to the peak power in the CIR
cir_mag2 = (abs(cir).^2);
cir_max = max(cir_mag2);
cir_mag2 = (abs(cir).^2)/cir_max;
nf = mean(cir_mag2(round(length(cir)*0.8):round(length(cir)*0.9)));
cir_peaks_k = find( cir_mag2 > max([10*nf cir_max/15.8489]) );
cir_peaks_t = t(cir_peaks_k);
cir_peaks = cir_mag2(cir_peaks_k);
if 1
plot(t, 10*log10(cir_mag2), cir_peaks_t, 10*log10(cir_peaks), 'ro')
refline(0,10*log10(nf)+10)
refline(0,-12)
end
if isempty(cir_peaks)
    K = NaN;
    return
end

% we want at least 10 components for computation.
% a component is considered 4 sample
if length(cir_peaks) < 40
    K = NaN;
    return
end

% calculate the peak value
[cir_max, cir_max_k] = max(cir_peaks);
cir_max_t = cir_peaks_t(cir_max_k);

% form the diffuse components to estimate nlos power
klos = cir_max_k;
dfuse_mag = cir_peaks; 
if klos > 1
    dfuse_mag(klos-1:klos+1) = NaN;
else
    dfuse_mag(klos:klos+1) = NaN;
end
dfuse_mag = dfuse_mag(~isnan(dfuse_mag));
dfuse_pwr = var(dfuse_mag);

% A. Doukas and G. Kalivas, "Rician K Factor Estimation for Wireless
% Communication Systems," 2006 International Conference on Wireless and
% Mobile Communications (ICWMC'06), Bucharest, 2006, pp. 69-69.  
% doi: 10.1109/ICWMC.2006.81  
K = 10*log10(cir_max^2/(2*dfuse_pwr));  

if K < 6
    K = NaN;
    return
end

% final filter: time-based determination.  The peak must occur within a
% reasonable amount of time from the start of the cir.  In this case, the
% calling function makes the time threshold determination.
cir_dT0 = cir_max_t - cir_peaks_t(1);
if cir_dT0 > dT0
    K = NaN;
    return
end
                
end

