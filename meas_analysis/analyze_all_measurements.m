% Analyze complex impulse responses from measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

clear;
clc;

%% cloud/stationary measurements
% cd('Boulder_c')
% delete('cloud_stats.mat')
% estimate_cloud_cwd('*.mat', [1 1 1 0 1]')
% cd('..')

%% mobile measurements
delete('stats.mat')
top_dirs = { 'AAplant', 'Boulder', 'GBurg'};
doall = true;
figvis = false;
for ii = 1:length(top_dirs)
    
    cd(top_dirs{ii})

    analyze_cwd('all');
    
    cd('..\')

end

