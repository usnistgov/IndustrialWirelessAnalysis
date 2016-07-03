function [ tau ] = compute_delay_spread( t, cir, minA )
% COMPUTE_DELAY_SPREAD Compute the delay spread of the input CIR
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

k = find(cir > minA);
cir_k = cir(k);
t_k = t(k);
t_k = t_k - t_k(1);  % remove prop delay
tau = t_k(:).'*abs(cir_k(:));  
tau = tau/sum(abs(cir_k));

end

