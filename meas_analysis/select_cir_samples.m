function [ k, nf ] = select_cir_samples( r, cir )
% SELECT_CIR_SAMPLES Compute the delay spread of the input CIR
%
% Outputs:
%   k is an array of the selected indices
%   nf is the noise floor
%
% Inputs:
%   r is the range in meters from transmitter to receiver
%   cir is the real or complex valued channel impulse response
%
% Time and CIR vectors must have the same length.
%
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

% only consider 
if r < 4 
    k = [];
    nf = nan;
    return
end

% linear domain thresholds
nft = 10;           % 10 dB above noise floor
pkt = 1/31.6228;    % 15 dB from peak

% eliminate the anomalous tail components (last 8 samples)
cir = cir(1:round(length(cir)*0.5));
% cir = cir(1:end-8);

% compute the magnitude squared of the cir normailzed to the peak value
cir_mag2 = (abs(cir).^2);
cir_max = max(cir_mag2);
cir_mag2 = (abs(cir).^2)/cir_max;

% select the cir sample for use in average cir estimation
nf = mean(cir_mag2(round(length(cir)*0.8):round(length(cir)*0.9)));
k = find( cir_mag2 > nf*nft & cir_mag2 > pkt );

end