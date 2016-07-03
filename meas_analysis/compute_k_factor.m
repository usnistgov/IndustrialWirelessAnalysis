function [ K ] = compute_k_factor( t, cir, r, L )
%COMPUTE_K_FACTOR Computes the K-factor of the impulse response
%   Computes the K factor of the CIR using all components of magnitude
%   greater than the noise floor times L.  The range, r, in meters guides
%   the extraction of the LOS power component.
%
%   t:      time in seconds
%   cir:    the channel impulse response
%   r:      range in meters 
%   L:      linear limit multiplier 
% 
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

% segregate the peak values
cir_mag = abs(cir);
nf = mean(cir_mag(round(length(cir)*0.8):round(length(cir)*0.9)));
cir_peaks = cir_mag(cir_mag > nf*L);
if isempty(cir_peaks)
    K = NaN;
    return
end

% find the peak based on the reported range from the measurements logs
% this will only work if the range estimate is accurate
C = 2.99792458e8; %m/s
Ts = t(2)-t(1);
klos = ceil((r/C)/Ts)+1;
los_peak = max(cir_mag(klos-1:klos+1));

% form the diffuse components to estimate nlos power
dfuse_mag = cir_peaks; 
if klos > 1
    dfuse_mag(klos-1:klos+1) = NaN;
else
    dfuse_mag(klos:klos+1) = NaN;
end
dfuse_pwr = nanvar(dfuse_mag);

% tpeaks = t(cir_mag>nf*L);
% cir_peaks = cir_mag(cir_mag>nf*L);
% plot(tpeaks*1e9, 10*log10(cir_peaks), tpeaks*1e9, 10*log10(dfuse_mag));
% drawnow

% A. Doukas and G. Kalivas, "Rician K Factor Estimation for Wireless
% Communication Systems," 2006 International Conference on Wireless and
% Mobile Communications (ICWMC'06), Bucharest, 2006, pp. 69-69.  
% doi: 10.1109/ICWMC.2006.81
K = 10*log10(los_peak^2/(2*dfuse_pwr));    
                
end

