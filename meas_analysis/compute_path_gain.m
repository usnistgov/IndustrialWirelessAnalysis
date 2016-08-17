function [ G ] = compute_path_gain( cir, Ga_tx, Ga_rx )
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

if isempty(cir)
    G = NaN;
    return
end

cir_mag2 = abs(cir).^2;
N = length(cir_mag2);
G = 10*log10(sum(cir_mag2)/N) - Ga_tx - Ga_rx; 

end

