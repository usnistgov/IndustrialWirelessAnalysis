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

% normalize the cir for computation
a_k = a_k/max(a_k);

% compute the peaks
findpeaks(a_k),xlim([1 45])
[pks, k_pks] = findpeaks(a_k);
t_k = t(k_pks);
t_k = t_k - t_k(1);  % remove propagation delay

% compute the mean excess delay
tau_u = sum(pks(:)'*t_k(:))/sum(pks.^2);

% compute second moment, rms delay spread
% tau_s = sqrt(sum(pks(:)'*t_k(:).^2)/sum(pks) - tau_u);
tau_s = sqrt( sum(pks(:)' * (t_k(:)-tau_u).^2) /  sum(pks)  );

% compute the duration
T = t_k(end) - t_k(1);

end

