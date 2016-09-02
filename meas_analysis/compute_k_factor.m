 function [ K, LOS, klos ] = compute_k_factor( t, cir, ns )
%COMPUTE_K_FACTOR Computes the K-factor of the impulse response
%   Computes the K factor of the CIR using all components of magnitude
%   greater than the noise floor times L.  The range, r, in meters guides
%   the extraction of the LOS power component.  It is assumed that the
%   select_cir_samples() function has been used to remove unwanted samples
%   withing the cir
%
%   Outputs:
%       K:      The K factor in dB
%       LOS:    boolean LOS indicator
%           -1 for NLOS
%            0 for unknown
%           +1 for LOS
%       k0 is the index of the first peak
%
%   Inputs:
%       t:      time in seconds
%       cir:    the channel impulse response
%       ns:    oversample rate
% 
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

K = NaN;
LOS = 0;
klos = nan;
if isempty(cir)
    return
end

% calculate the peak value and the time of its occurrance
cir_mag = abs(cir);
cir_max = max(cir_mag);

% Time-based determination.  The first peak must be greater in magnitude
% than all of the other peaks in the cir.
[pks, k_pks] = findpeaks(cir_mag,1:length(cir_mag));  % requires the signal processing toolbox
if isempty(pks)
    return
end
klos = k_pks(1);
cir_los = cir_mag(klos);
if cir_los ~= cir_max
    LOS = -1;  % max is not the first peak, so NLOS
    return;
end

% form the diffuse components to estimate non-los power.  This is done by
% removing the local samples related to peak magnitude value
dfuse_mag = cir_mag; 
if klos > ns/2
    kklos = klos-1:klos+1;
    dfuse_mag(kklos) = NaN;
else
    kklos = klos:klos+1;
    dfuse_mag(kklos) = NaN;
end
dfuse_mag = dfuse_mag(~isnan(dfuse_mag));
dfuse_pwr = var(dfuse_mag);

% A. Doukas and G. Kalivas, "Rician K Factor Estimation for Wireless
% Communication Systems," 2006 International Conference on Wireless and
% Mobile Communications (ICWMC'06), Bucharest, 2006, pp. 69-69.  
% doi: 10.1109/ICWMC.2006.81  
LOS = 1;
K = 10*log10(cir_los^2/(2*dfuse_pwr));  
                
end

