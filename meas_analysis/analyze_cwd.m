% Analyze complex impulse responses from measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

top_dirs = { 'AAPlant', 'Boulder', 'GBurg'};
doall = true;
figvis = false;

for ii = 1:length(top_dirs)
    
    cd(top_dirs{ii})

    set(0,'DefaultFigureWindowStyle','docked')
    %set(0,'DefaultFigureWindowStyle','normal')

    pattern = '*.mat';
    if exist('TEST_DATA','var')
        estimate_channel_cwd(pattern, doall, TEST_DATA);
    else
        estimate_channel_cwd(pattern, doall);
    end
    
    cd('..\')

end

