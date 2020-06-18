function [ k, nf, cir1, pk_pwr ] = select_cir_samples( r, cir )
% SELECT_CIR_SAMPLES Compute the delay spread of the input CIR
%
% Outputs:
%   k is an array of the selected indices
%   nf is the noise floor
%   cir1 is the modified cir, noise set to zero
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

k = [];
nf = nan;
cir1 = nan;
pk_pwr = nan;

% % only consider non-near measurements
% if r < 3
%     return;
% end

% estimate peak power
pk_pwr = 10*log10(abs(max(cir).^2));

% linear domain thresholds
nft = 10;           % 10 dB above noise floor
pkt = 1/1000.0;    % 10-20 dB from peak is ITU recommendation
% the dynamic range of our instrumentation allows for a larger calculation
% domain but we select 30 dB here but the delay spread calculation we
% further refine to 15 dB

% eliminate the anomalous tail components 
cir = cir(1:round(length(cir)*0.8));

% compute the magnitude squared of the cir normailzed to the peak value
cir_mag2 = (abs(cir).^2);
cir_max = max(cir_mag2);
cir_mag2 = (abs(cir).^2)/cir_max;

% select the cir sample for use in average cir estimation
nf = mean(cir_mag2(round(length(cir)*0.8):round(length(cir)*1.0)));
k = find( (cir_mag2 > nf*nft) & (cir_mag2 > pkt) );

% modified cir
cir1 = zeros(size(cir));
cir1(k) = cir(k);

%plot([10*log10(cir_mag2) 10*log10(abs(cir1/max(abs(cir1))).^2)],'+-'),xlim([0 500]),refline(0,10*log10(nf)+10),refline(0,-30)

end

