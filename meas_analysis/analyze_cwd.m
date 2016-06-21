% Analyze complex impulse responses from measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

addpath('C:\git\wpsm_rf_analysis')

tic
doall = false;
figvis = false;

set(0,'DefaultFigureWindowStyle','docked')
%set(0,'DefaultFigureWindowStyle','normal')

pattern = '*.mat';
estimate_channel_cwd(pattern, doall, figvis, 12.5E-9);

toc