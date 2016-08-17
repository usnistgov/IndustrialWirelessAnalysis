function [ K, LOS ] = compute_k_factor( t, cir, ns )
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
%
%   Inputs:
%       t:      time in seconds
%       cir:    the channel impulse response
%       ns:    oversample rate
% 
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

if isempty(cir)
    K = NaN;
    LOS = 0;
    return
end

% we want at least 4 components for computation.
% a component is considered to be at least 3 samples
% Note that this is a rought estimate of the number of components for K
% factor estimation.  A better method could be to use the delay spread to
% determine candidacy for K estimation.
if length(cir) < 12
    K = NaN;
    LOS = 0;
    return
end

% default is LOS
LOS = +1;

% calculate the peak value
cir_mag = abs(cir);
[cir_max, cir_max_k] = max(cir_mag);
cir_max_t = t(cir_max_k);

% form the diffuse components to estimate non-los power.  This is done by
% removing the local samples related to peak magnitude value
klos = cir_max_k;
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
K = 10*log10(cir_max^2/(2*dfuse_pwr));  

if K < 6
    K = NaN;
    LOS = -1;
    return
end

% final filter: time-based determination.  The peak must occur within a
% reasonable amount of time from the start of the cir.  In this case, the
% calling function makes the time threshold determination. 
Ts = t(2)-t(1);
dT0 = ns*Ts;
cir_dT0 = cir_max_t - t(1);
if cir_dT0 > dT0
    K = NaN;
    LOS = -1;
    return
end
                
end

