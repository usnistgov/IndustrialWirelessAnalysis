function [ c_k ] = select_for_avg_cir( cir )
% SELECT_FOR_AVG_CIR Select CIR components for average cir computation.
% Propagation delay is removed.
%
% Outputs:
%   c_k is the continuous time-normalized cir (complex)
%
% Inputs:
%   cir is the real or complex valued channel impulse response
%
% Time and CIR vectors must have the same length.
%
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

if isempty(cir)
    c_k = 0;
    return
end

% we want at least 4 components for computation.
% a component is considered to be at least 3 samples
% Note that this is a rought estimate of the number of components for K
% factor estimation.  A better method could be to use the delay spread to
% determine candidacy for K estimation.
if length(cir) < 12
    c_k = 0;
    return
end

% normalize the cir for computation
c_k = cir/max(abs(cir));

% remove last half of cir
c_k = c_k(1:round(0.5*(length(c_k))));

% compute the noise floor of the sample
nf_domain = round(length(c_k)/4);
nf = sum(abs(c_k(end-nf_domain:end)).^2)/nf_domain;

% only consider cir's with GT 20 dB SNR peak to noise floor
if 10*log10(max(abs(c_k).^2)/nf) < 30
    c_k = 0;
    return
end

% Zero out the noise floor and remove the samples before the beginning of the
% CIR. This removes the propagation delay but preserves the delay spread
% information.  We only consider samples 10 dB above the noise floor.  The
% Peak-Noise SNR is already determined to be at least 30 dB, so the number
% of samples will be non-zero.
nf_thresh = 10;
k_zero = abs(c_k).^2 < nf*nf_thresh;
c_k( k_zero ) = 0;
first_k = find(c_k,1,'first');
if first_k > 1
    c_k(1:first_k-1) = [];
else
    c_k(1) = [];
end

end

