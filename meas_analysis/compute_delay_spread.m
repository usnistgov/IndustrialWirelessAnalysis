function [ tau_u, tau_s, T ] = compute_delay_spread( t, cir, nf )
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
%   nf is the linear domain measured noise floor
%
% Time and CIR vectors must have the same length.
%
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

if isempty(cir)
    tau_u = NaN;
    tau_s = NaN;
    T = NaN;
    return
end

% we want at least 4 components for computation.
% a component is considered to be at least 3 samples
% Note that this is a rought estimate of the number of components for K
% factor estimation.  A better method could be to use the delay spread to
% determine candidacy for K estimation.
if length(cir) < 12
    tau_u = NaN;
    tau_s = NaN;
    T = NaN;
    return
end

a_k = abs(cir);

% only consider cir's with GT 20 dB SNR peak to noise floor
% if 10*log10(max(a_k.^2)/nf) < 20
%     tau_u = NaN;
%     tau_s = NaN;
%     T = NaN;
%     return
% end

% normalize the cir for computation
a_k = a_k/max(a_k);
t_k = t;
t_k = t_k - t_k(1);  % remove propagation delay

% remove spurious outliers
t_k_x = 3*median(t_k);
a_k = a_k(t_k<t_k_x);
t_k = t_k(t_k<t_k_x);

% compute the power in the signal
% P_k = abs(a_k).^2;
P_k = a_k.^2;

% compute the mean excess delay
tau_u = sum(P_k(:)'*t_k(:))/sum(P_k);  

% compute second moment, rms delay spread
tau_s = sqrt(sum(P_k(:)'*t_k(:).^2)/sum(P_k) - tau_u^2);

% compute the duration
T = t_k(end) - t_k(1);

end

