% Convert files to struct of meta + array of CIR's
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

clear;
clc;
%addpath('C:\git\wpsm_rf_analysis')
dirs = { 'Gaithersburg Day 1', 'Gaithersburg Day 2', 'Gaithersburg Day 3', 'AAPlant D1', 'AAPlant D2', 'AAPlant D3' };
pattern = '*.mat';

dbstop error
for kk = 1:length(dirs)
    disp(['entering ' dirs{kk}])
    chdir(dirs{kk})
    convert_to_array(pattern);
    chdir('..')
    disp(['now in ' chdir()])
end




