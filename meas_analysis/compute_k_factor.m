function [ K ] = compute_k_factor( t, cir, r, L, dT0 )
%COMPUTE_K_FACTOR Computes the K-factor of the impulse response
%   Computes the K factor of the CIR using all components of magnitude
%   greater than the noise floor times L.  The range, r, in meters guides
%   the extraction of the LOS power component.
%
%   t:      time in seconds
%   cir:    the channel impulse response
%   r:      range in meters (not used)
%   L:      linear limit multiplier 
%   dT0:    time to expect peak from first cir sample
% 
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

C = 2.99792458e8; %m/s
Ts = t(2)-t(1);

% eliminate the anomolous tail components (last 8 samples)
t = t(1:end-8);
cir = cir(1:end-8);

% segregate the peak values above the noise floor and within 10 dB of the
% peak magnitudea
cir_mag = abs(cir);
nf = mean(cir_mag(round(length(cir)*0.8):round(length(cir)*0.9)));
cir_max = max(cir_mag);
cir_peaks_k = find( cir_mag > max([nf*L cir_max/10]) );
cir_peaks_t = t(cir_peaks_k);
cir_peaks = cir_mag(cir_peaks_k);
if isempty(cir_peaks)
    K = NaN;
    return
end

% we want at least 10 components for computation.
if length(cir_peaks) < 10
    K = NaN;
    return
end

%plot(cir_peaks_t, 10*log10(cir_peaks),'x')

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
dfuse_pwr = nanvar(dfuse_mag);

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

