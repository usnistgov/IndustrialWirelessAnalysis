function [ tau_u, tau_s, T ] = compute_delay_spread(obj, Ts, cir )
% COMPUTE_DELAY_SPREAD Compute the delay spread of the input CIR
%
% Outputs:
%   tau_u is the mean excess delay
%   tau_s is the rms excess delay
%   T is the delay interval
%
% Inputs:
%   Ts is the sample interval of the cir
%   cir is the real or complex valued channel impulse response
%
% Time and CIR vectors must have the same length.
%
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

tau_u = NaN;
tau_s = NaN;
T = NaN;

if isempty(cir)
    return
end

% we want at least 4 components for computation.
% a component is considered to be at least 3 samples
% Note that this is a rought estimate of the number of components for K
% factor estimation.  A better method could be to use the delay spread to
% determine candidacy for K estimation.
% if length(cir) < 12
%     return
% end

t = Ts*(0:length(cir)-1);
t = t - t(1);  

% compute the cir magnitude
cir_mag = abs(cir);
p = cir_mag.^2;
p = p/max(p);

% delay spread sensitivity to spurious outliers
pkt = 1/31.6;    % 15 dB from peak is ITU-R 1407-5 recommendation
p(p<pkt) = 0.0;

% compute the mean excess delay
tau_u = sum(p(:)'*t(:))/sum(p);

% compute second moment, rms delay spread
% tau_s = sqrt(sum(pks(:)'*t_k(:).^2)/sum(pks) - tau_u);
tau_s = sqrt( sum(p(:)' * (t(:)-tau_u).^2) /  sum(p)  );

% compute the delay interval
k_t = find(cir_mag>0);
T = t(k_t(end));

end

