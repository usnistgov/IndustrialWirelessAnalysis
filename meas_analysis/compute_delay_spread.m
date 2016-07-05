function [ tau_u, tau_s, T ] = compute_delay_spread( t, cir, minA )
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
%   minA is the minimum amplitide to be considered for computation with
%   units of linear amplitude
%
% Time and CIR vectors must have the same length.
%
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

if nargin < 3
    minA = -Inf;
end

% select components above noise threshold
k = find(cir > minA);

% eliminate the last 10% of the CIR do to fft wrapping from post-processing
% note that this wrapping is not ideal, but a result of the method of
% post-processing of the raw data.
NN = size(cir);
NN = round(NN*0.9);

% select components
a_k = cir(k(1:NN));
t_k = t(k(1:NN));
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

