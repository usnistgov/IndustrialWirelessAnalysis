% Analyze complex impulse responses from measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

addpath('C:\git\wpsm_rf_analysis')

tic
doall = true;
figvis = false;

set(0,'DefaultFigureWindowStyle','docked')
%set(0,'DefaultFigureWindowStyle','normal')

pattern = '*.mat';
if exist('TEST_DATA','var')
    estimate_channel_cwd(pattern, doall, figvis, TEST_DATA);
else
    estimate_channel_cwd(pattern, doall, figvis);
end

toc