function [ k ] = select_cir_samples( r, t, cir )
% SELECT_CIR_SAMPLES Compute the delay spread of the input CIR
%
% Outputs:
%   k is an array of the selected indices
%
% Inputs:
%   r is the range in meters from transmitter to receiver
%   t is the time vector
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
    return
end

% linear domain thresholds
nft = 10;           % 10 dB above noise floor
pkt = 1/31.6228; 1/15.8489;    % 12dB from peak

% eliminate the anomalous tail components (last 8 samples)
cir = cir(1:end-8);

% compute the magnitude squared of the cir normailzed to the peak value
cir_mag2 = (abs(cir).^2);
cir_max = max(cir_mag2);
cir_mag2 = (abs(cir).^2)/cir_max;
nf = mean(cir_mag2(round(length(cir)*0.8):round(length(cir)*0.9)));
k = find( cir_mag2 > max([nft*nf pkt]) );

if 0
    t = t(1:end-8);
    cir_peaks = cir_mag2(k);
    cir_peaks_t = t(k);
    plot(t, 10*log10(cir_mag2), cir_peaks_t, 10*log10(cir_peaks), 'ro')
    refline(0,10*log10(nft*nf))
    refline(0,10*log10(pkt))
end


end