% IQdata_CloudLocations and Strct_Metadata should be loaded from the *.mat
% files for each of the smart manufacturing measurements. Run this after
% loading the *.mat file.

% After this script is run, the positions for each record will be held in 
% the 'xpos' and 'ypos' variables. For example, IQdata(:,1) corresponds to
% xpos(1) and ypos(1).

IQdata_CloudLocations_m_num = DATA.IQdata_CloudLocations_m_num;
Strct_Metadata = DATA.Strct_Metadata;

% Get x Positions and repeat each element by the number of records per acquisition
xpos = repelem(IQdata_CloudLocations_m_num.xPositions,Strct_Metadata.NumberRecordperAcqusition_num,1);
% Get y Positions and repeat each element by the number of records per acquisition
ypos = repelem(IQdata_CloudLocations_m_num.yPositions,Strct_Metadata.NumberRecordperAcqusition_num,1);


%plot the positions
%{
fig = figure();
scatter(xpos,ypos);
xlabel('X positions (m)'); ylabel('Y positions (m)');
%}