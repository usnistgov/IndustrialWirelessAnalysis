% Analyze complex impulse responses from measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

clear;
clc;

top_dirs = { 'AAPlant', 'Boulder', 'GBurg'};
doall = true;
figvis = false;

for ii = 1:length(top_dirs)
    
    cd(top_dirs{ii})

    analyze_cwd
    
    cd('..\')

end

