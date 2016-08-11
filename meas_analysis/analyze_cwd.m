% Analyze complex impulse responses from measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

doall = true;
figvis = false;

set(0,'DefaultFigureWindowStyle','docked')
%set(0,'DefaultFigureWindowStyle','normal')

pattern = '*.mat';
if exist('TEST_DATA','var')
    estimate_channel_cwd(pattern, doall, TEST_DATA);
else
    estimate_channel_cwd(pattern, doall);
end


