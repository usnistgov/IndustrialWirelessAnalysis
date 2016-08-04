function [ tau_u, tau_s, T ] = compute_delay_spread( t, cir )
% COMPUTE_DELAY_SPREAD Compute the delay spread of the input CIR
%
% Outputs:
%   tau_u is the mean excess delay
%   tau_s is the rms excess delay
%   T is the duration of the CIR
%
% Inputs:
%   t is the time vector
%   cir is the real or complex valued channel impulse response
%
% Time and CIR vectors must have the same length.
%
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

if nargin < 3
    minA = -Inf;
end

% compute the magnitude of the cir
cir_mag = abs(cir);

% select components above noise threshold and within 10 dB of the peak
cir_max = max(cir_mag);
k = find(cir_mag > max([minA cir_max/10]));

% select components
a_k = cir_mag(k);
t_k = t(k);
t_k = t_k - t_k(1);  % remove prop delay

% compute the power in the signal
P_k = abs(a_k).^2;

% compute the mean excess delay
tau_u = sum(P_k(:)'*t_k(:))/sum(P_k);  

% compute second moment, rms delay spread
tau_s = sqrt(sum(P_k(:)'*t_k(:).^2)/sum(P_k) - tau_u^2);

% compute the duration
T = t_k(end) - t_k(1);

end

