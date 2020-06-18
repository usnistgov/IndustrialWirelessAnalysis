function [ g_dB, g ] = compute_path_gain( cir, Ga_tx, Ga_rx )
% COMPUTE_PATH_GAIN Compute the path gain of the transmission
%
% Outputs:
%   G is the power gain in dB
%   Ga_tx, Ga_rx are the antenna gains in dBi
%
% Inputs:
%   cir is the real or complex valued channel impulse response
%
% Time and CIR vectors must have the same length.
%
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

g_dB = NaN;
g = NaN;
if isempty(cir)
    return
end

% cir_mag2 = cir(:)'*cir(:);
% g_dB = 10*log10(cir_mag2) - Ga_tx - Ga_rx;
% g = power(10, g_dB/10);

MM = length(cir);
MMd2 = floor(MM/2);
th = 1.0000e-11; %10^(-110/10);

% compute mag squared
mag2 = abs(cir).^2;

% remove wrapping if it exists
Lsh = MMd2+find(mag2(MM/2:end)>th)-1;
if ~isempty(Lsh) % remove wrapping effect
    Lsh = MM-Lsh(1);
    cir = circshift(cir,Lsh);
    mag2 = abs(cir).^2;
end
mag2(mag2<th) = [];  % SNR validation

g_dB = 10*log10(sum(mag2)) - Ga_tx - Ga_rx;
g = power(10, g_dB/10);

end

