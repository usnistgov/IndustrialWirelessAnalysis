 function [ K, LOS, k_pks ] = compute_k_factor( cir, ns )
%COMPUTE_K_FACTOR Computes the K-factor of the impulse response
%   Computes the K factor of the CIR using all components of magnitude
%   greater than the noise floor times 10.
%
%   Outputs:
%       K:      The K factor in dB
%       LOS:    boolean LOS indicator
%           -1 for NLOS
%            0 for unknown
%           +1 for LOS
%       k_pks are the peak indices
%
%   Inputs:
%       cir:    the channel impulse response
%       ns:    oversample rate
% 
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

K = NaN;
LOS = 0;
k_pks = nan;
if isempty(cir)
    return
end

% normalize the cir
cir = cir/max(abs(cir));

% calculate the peak value and the time of its occurrance
cir_mag2 = abs(cir).^2;
cir_max2 = max(cir_mag2);

% Time-based determination.  The first peak must be greater in magnitude
% than all of the other peaks in the cir.
[pks, k_pks] = findpeaks(cir_mag2,1:length(cir_mag2)); 
if isempty(pks)
    return
end

% the los peak must be at or after the PN oversample factor
[~, max_i] = max(pks);
if max_i > 1
    LOS = -1;  % max is not the first peak, so NLOS
    return;    
end
% kpeak = k_pks(1);
% cir_los = cir_mag2(kpeak);
% if cir_los ~= cir_max2
%     LOS = -1;  % max is not the first peak, so NLOS
%     return;
% end

% JPL's Wireless Communication Reference Website
% Chapter: Wireless Channels
% Section Rician Channels, Indoor Channels
% URL http://www.wirelesscommunication.nl/reference/chaptr03/ricepdf/measrice.htm
K = 10*log10(pks(1)/mean(pks(2:end)));
LOS = +1;
                
end

