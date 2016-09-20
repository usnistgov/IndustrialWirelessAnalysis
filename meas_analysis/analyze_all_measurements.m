% Analyze complex impulse responses from measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

clear;
clc;

%% cloud/stationary measurements
cd('Boulder_c')
estimate_cloud_cwd('*.mat', ones(3,1))
cd('..')

%% mobile measurements
top_dirs = { 'AAplant', 'Boulder', 'GBurg'};
doall = true;
figvis = false;
for ii = 1:length(top_dirs)
    
    cd(top_dirs{ii})

    analyze_cwd('all');
    
    cd('..\')

end

