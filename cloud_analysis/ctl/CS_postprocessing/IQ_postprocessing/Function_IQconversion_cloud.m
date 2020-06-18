%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%This function saves the raw .tdms file into a .mat file.  It also
%converts the .tdms file into Smeas/Sref*attenuation
%Assumptions:
% Function: Function_IQconversion
% Author: Jeanne Quimby - 3/17/2016
%         David Novotny
%         Alexandra Curtin
%--------------------------------------------------------------------------
%Inputs and Outputs for Function_IQconversion
%Inputs:
%   * none
%Outputs
%   * text file
%--------------------------------------------------------------------------
%Function calls
%    1) convertTDMS - converts labview TDMS files to matlab cells
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function Function_IQconversion_cloud(MeasurementPath,t_file, t_arr, t_rec, Collection, header, ...
    numFiles, Tx_wfm, normalize,filter_IQdata,Strct_Metadata, IQdata_CloudLocations, min_amplitude)

MatFile_str = Strct_Metadata.MatFile_str;                   %Save Filename
str_savefile = [MeasurementPath MatFile_str '_pp.mat'];
str_savefileforus = [MeasurementPath MatFile_str '_ppforNIST.mat'];
if (exist(str_savefile) ~= 0)
       delete(str_savefile);
end
    
Strct_Metadata.IQdata_Timing_Header_str = ['Acqusition Index' 'Record Index' 'Acqusition Time(date and seconds)' 'Record Time (seconds)' 'Total Time (date and seconds)'];
Strct_Metadata.IQdata_Range_Header_str = ['Acqusition Index' 'Record Index' 'Range Distance(m)'];
Strct_Metadata.IQdata_Header_str = 'IQ Data (real and complex';

numR=length(t_arr);  wordL=length(t_rec);

%modifying to avoid any possibility of mistake
indexertemp = 0;


%My modifications - Jdiener

%Create Frequency Conversion
if isfield(Strct_Metadata, 'Frequency_GHz_num')
    freqc = Strct_Metadata.Frequency_GHz_num;
    freqc = freqc*1e9;
    dt = t_rec(2) - t_rec(1); df = 1/dt;
    freqvec = transpose(-df/2 : df/size(filter_IQdata, 1) : df/2 - 1) + freqc;
    freqvec = freqvec./1e9; %convert to units of GHz.
    freqvector = freqvec.*1e9;
    
end

%%%Pre index the IQ data set - this should be slightly faster. 
raw_IQdata = zeros(Strct_Metadata.CodewordLength_num, ...
    Strct_Metadata.NumberAcqusitions_num*Strct_Metadata.NumberRecordperAcqusition_num );

IQdata_scaled = zeros(Strct_Metadata.CodewordLength_num, ...
    Strct_Metadata.NumberAcqusitions_num*Strct_Metadata.NumberRecordperAcqusition_num );

% % % % % pathlossfreqrecavg = zeros(length(freqvec), 1); %%%This takes the Frequency data from every record, constructs a runniung average,
% % % % % %%%and uses this to estimate the path loss over FREQUENCY. Note that the
% % % % % %%%averaging is performed using the magnitude, NOT the complex, frequency
% % % % % %%%response. This is NOT square filtered, and NOT Filtered by the PN
% % % % % %%%filter. 

%%%CHANGING TO - 
FreqResponse_allfreq = zeros(Strct_Metadata.CodewordLength_num, ...
    Strct_Metadata.NumberAcqusitions_num*Strct_Metadata.NumberRecordperAcqusition_num );



PathGain_avgfreq = zeros(Strct_Metadata.NumberAcqusitions_num*Strct_Metadata.NumberRecordperAcqusition_num, 1);
%%%This contains the pathloss Frequency estimate. What this does is take
%%%the frequency data, filter it, and then take the power of this data
%%%sum(abs(data)). This is then scaled by a factor which undoes the power
%%%of the filter in an average sense (i.e. if the filter on average removes
%%%power, the scaling factor will add power back in). 

PathGain_avgtime = PathGain_avgfreq; 
%%%This contains the Pathloss Area ESTIMATE, as obtained from
%%%%analyze record. This is done by taking the sum of the magnitude of the
%%%%CIR data. This uses the square_filtered frequency domain data, and will
%%%%represent a slightly different estimate from the time domain PDP path
%%%%loss calculation methods.


%%%Vector Grab / Frequency Limited (Through Square-Filter)
%%%Concern ourselves with only those frequency points that have relative
%%%energy above some value, as judged by the PN spectrum. Thus, if we grab
%%%points with realtive amplitude above 0.1, then...
% % % min_amplitude = 0.1; - Grabbed from Input. 0.1 is a good value.
vectorgrab = find(filter_IQdata >= min_amplitude);
%%%%

indexer = 0; %Linear index for writing data into matrices. More convenient than the original representation.



for j=1:numFiles
    
    for c=1:Collection{j,4}{11,2}
        
        for t=1:numR
            
            indexer = indexer + 1;
            
            dummyRec=getRecord(Collection{j,1},c,t,wordL,numR);         
            
            [IQ_data, CIR_Discretescaling_Bandlimit, FreqResponse_allfreqtmp, PathGain_avgfreqtmp, PathGain_avgtimetmp]= ...
            analyzeRecord_vs5(dummyRec,Tx_wfm,t_rec,normalize,freqvector,[],filter_IQdata, vectorgrab);        
        
            raw_IQdata(:,indexer)= IQ_data;
            IQdata_scaled(:, indexer) = CIR_Discretescaling_Bandlimit;
            FreqResponse_allfreq(:,indexer) = FreqResponse_allfreqtmp;            
            PathGain_avgfreq(indexer) = PathGain_avgfreqtmp;
            PathGain_avgtime(indexer) = PathGain_avgtimetmp;
            
            
            
        end
    end
end


%%%%Compute RMS Delays

[~, rmsvalues] = RMSdelay( abs(IQdata_scaled).^2 ,t_rec, 30 ) ;  %threshhold in db - timevector T

figure; plot(rmsvalues);
title('RMS Values');
ylabel('RMS Spread (ns)');
grid on;

% keyboard;

%%%%

pathlossavgrec_dB =  10.*log10(PathGain_avgtime) - Strct_Metadata.ReceiverAntennaGain_dBi_num ...
                                 - Strct_Metadata.TransmitterAntennaGain_dBi_num;

%Compute Mean and Std
for j = 1:numFiles;
    ndxj = 1+(j-1).*numR:numR.*j;
    meanPathGain(j) = sum(PathGain_avgtime(ndxj))./numR;
    varPathGain(j) = sum((PathGain_avgtime(ndxj) - meanPathGain(j)).^2)./(numR - 1);
end    
stdPathGain = sqrt(varPathGain);
stdPathGain_dB = 10.*log10(abs(stdPathGain+meanPathGain)./meanPathGain);    
meanPathGain_dB =  10.*log10(meanPathGain) - Strct_Metadata.ReceiverAntennaGain_dBi_num ...
                                 - Strct_Metadata.TransmitterAntennaGain_dBi_num;

 %%%Compute MEAN RMS
for j = 1:numFiles;
    ndxj = 1+(j-1).*numR:numR.*j;
    meanRMS(j) = sum(rmsvalues(ndxj))./numR;
    varrms(j) = sum((rmsvalues(ndxj) - meanRMS(j)).^2)./(numR - 1);
end                                
 
                             
%meanPathGain_dB = [ones(1,20).*5 ones(1,20).*10  1:9]
IQdata_Timing_cll = (1:size(IQ_data,1))';

% Compute distance
IQdata_Range_m(:,6) = Strct_Metadata.Tx_xyz_m_cll{3} - Strct_Metadata.Rx_xyz_m_cll{3};
[XPosition,YPosition,metersDistance,meanPathGain_dBp,meanPathGain_dBl] = ...
    distanceFromCenterOfCloud(Strct_Metadata,Strct_Metadata.Tx_xyz_m_cll,Strct_Metadata.Rx_xyz_m_cll,IQdata_CloudLocations,meanPathGain);

% keyboard

[~,~,~,rmsdata_spatial] = distanceFromCenterOfCloud_RMS((Strct_Metadata.Tx_xyz_m_cll), ...
    (Strct_Metadata.Rx_xyz_m_cll), IQdata_CloudLocations,meanRMS);

% Plot the actual data.
figure
errorbar(metersDistance, ... % X-axis
    meanPathGain_dB,stdPathGain_dB, ... % Upper error bar
    'Marker', 's', 'MarkerSize', 3, ...
    'LineStyle', 'none','LineWidth',2);
str_title =MatFile_str; ndx = strfind(str_title,'_'); str_title(ndx) = ' ';
title(str_title); set(gca,'FontSize',30); set(gca,'FontName','times')      
xlabel('Range (meters)');  ylabel('Path Gain (dB)');

meanPG_values = unique(round(meanPathGain_dB));

%%%mean rms
meanrms_values = unique(round(meanRMS));

%Plot 3D
%figure
%surface(XPosition,YPosition,meanPathGain_dBp,'EdgeColor','none'); colorbar('Limits',[min(meanPG_values) max(meanPG_values)],'Ticks',meanPG_values);
%hold on
%plot(Strct_Metadata.Rx_xyz_m_cll{1} - IQdata_CloudLocations(:,1),Strct_Metadata.Rx_xyz_m_cll{2} - IQdata_CloudLocations(:,2),'k*')
%title([str_title ' spline']); set(gca,'FontSize',30); set(gca,'FontName','times');
%xlabel('x position (meters)'); ylabel('y position (meters)');  zlabel('Path Gain (dB)');

%figure; surface(XPosition,YPosition,meanPathGain_dBl,'EdgeColor','none'); colorbar('Limits',[min(meanPG_values) max(meanPG_values)],'Ticks',meanPG_values);
%hold on
%plot(Strct_Metadata.Rx_xyz_m_cll{1} - IQdata_CloudLocations(:,1),Strct_Metadata.Rx_xyz_m_cll{2} - IQdata_CloudLocations(:,2),'k*')
%title([str_title ' linear']); set(gca,'FontSize',30); set(gca,'FontName','times');
%xlabel('x position (meters)'); ylabel('y position (meters)');  zlabel('Path Gain (dB)');

%figure; contourf(XPosition,YPosition,meanPathGain_dBp); colorbar('Limits',[min(meanPG_values) max(meanPG_values)],'Ticks',meanPG_values);
%hold on
%plot(Strct_Metadata.Rx_xyz_m_cll{1} - IQdata_CloudLocations(:,1),Strct_Metadata.Rx_xyz_m_cll{2} - IQdata_CloudLocations(:,2),'k*')
%title([str_title ' spline']); set(gca,'FontSize',30); set(gca,'FontName','times')      
%xlabel('x position (meters)'); ylabel('y position (meters)');  zlabel('Path Gain (dB)');

figure;
contourf(XPosition,YPosition,meanPathGain_dBl);
%contour3(YPosition,XPosition,meanPathGain_dBl,15,'k'); hold on; surf(YPosition,XPosition,meanPathGain_dBl, 'Edgecolor', 'none');
h = colorbar('Limits',[min(meanPG_values) max(meanPG_values)],'Ticks',meanPG_values); ylabel(h,'Path Gain (dB)','Rotation',270)
hold on
plot(Strct_Metadata.Rx_xyz_m_cll{1} - IQdata_CloudLocations(:,1),Strct_Metadata.Rx_xyz_m_cll{2} - IQdata_CloudLocations(:,2),'k*')
title([str_title ' linear']); set(gca,'FontSize',30); set(gca,'FontName','times')      
xlabel('x position (meters)'); ylabel('y position (meters)');  zlabel('Path Gain (dB)');

folderpath = regexprep(str_title,' ','','emptymatch') ;

%%%%Make a sub folder for each of the data runs - 
%jpegsavepath = '\\cfs2w.nist.gov\672\public\Quimby\Students\JDiener\Cloud_RMS_Temp\';
jpegsavepath = './figs'
jpegsavepathbase = fullfile(jpegsavepath,folderpath);

folderexist = exist(jpegsavepathbase); %%%Check if the folder exists - if 0, it does NOT
if folderexist == 0
    mkdir(jpegsavepathbase); %%%IF zero, MAKE folder. 
end
%%%%%




%%%%save the figures for the RMS Delay spread...





%set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.


% sefsef

%jpegsavestring = savestring(1:end - 4);
% % % jpegsavepath = '\\cfs2w.nist.gov\672\public\Quimby\Students\JDiener\Cloud_RMS_Temp\';

figure;
contourf(XPosition,YPosition,rmsdata_spatial);
hold on
plot(Strct_Metadata.Rx_xyz_m_cll{1} - IQdata_CloudLocations(:,1),Strct_Metadata.Rx_xyz_m_cll{2} - IQdata_CloudLocations(:,2),'k*')
xlabel('x position (meters)'); ylabel('y position (meters)');  zlabel('RMS Delay (ns)');
title({['RMS Delay Spread'], [str_title]});
% xlim([min(XPosition)-0.0001 max(XPosition)+0.0001])
% 
% % % %  ylim([min(YPosition)-0.001 max(YPosition)+0.0001])
% % % %  h = gca;
% % % %  qq = h.YTickLabel;
% % % %  dtick = (str2num(qq{2})) - str2num(qq{1});
% % % %  minval = str2num(qq{1}); maxval = str2num(qq{end});
% % % %  newmin = minval - dtick; newmax = maxval + dtick;
% % % %  
% % % %  
% % % %   ww = h.YTick;
% % % %  ww = [newmin  ww  newmax];
% % % %  
% % % %  set(h, 'YTick', ww)
% % % %  %%%now...
% % % %  ww = [num2str(newmin) ; qq ; num2str(newmax)];
% % % %  set(h, 'YTickLabel', ww);
% % % %  

 

%%%For ticks, do no more than 10...
%%%fit to units of 5....
fitvector = 20:5:125;
%%%now, do - 
[~, maxrmsindex] = min( abs( fitvector - (max(meanrms_values) + 0) )  ); 
%%%then grab the next highest one....
[~, minrmsindex] = min( abs( fitvector - (min(meanrms_values) - 0) )  ); 
if maxrmsindex < length(fitvector)
    topgrab = maxrmsindex + 1;
else
    topgrab = maxrmsindex;
end
if minrmsindex == 1
    botgrab = 1;
else
    botgrab = 2;
end
    
    
dtick = (fitvector(topgrab) - fitvector(botgrab))./ 10;
%%%rmstickvector 
rmstickvector = ((1:11)-1).*dtick + fitvector(botgrab) ;


%h = colorbar('Limits',[min(meanrms_values)-0 max(meanrms_values) + 0],'Ticks',meanrms_values);

%set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.

h = colorbar('Limits',[fitvector(botgrab) fitvector(topgrab)],'Ticks', rmstickvector);

 ylabel(h,'RMS Delay (ns)','Rotation',270)
%ylabel(h,'RMS Delay (ns)' );
 %ylabh = get(h,'YLabel');
 ylabh = get(h, 'Label');
 pause(0.0000001)
 pos1 = get(ylabh, 'Position');
 
 
 set(ylabh, 'Position', get(ylabh,'Position') + [2.6 0 0 ] )
pos2 = get(ylabh, 'Position');
%  set(ylabh,'Position',get(ylabh,'Position') + [1.5 0.0 0])



set(gca, 'FontSize', 24)

set(gcf, 'PaperUnits', 'normalized');
set(gcf, 'PaperPosition', [0 0 2 1]); %x_width=10cm y_width=15cm
% saveas(gcf, jpegsavepath, 'jpeg');


%jpegsavepath2 = [jpegsavepath 'new.jpg'];
jpegsavestring = ['RMSDelayPlotsContour'];
jpegsavepath = [jpegsavepathbase '\' jpegsavestring];
jpegsavepath2 = jpegsavepath;


%%%Apparently a more exact way.
pause(0.0001)
hgexport(gcf, jpegsavepath2,  ...
    hgexport('factorystyle'), 'Format', 'jpeg');

% keyboard

%%%%Another RMS Plot...

% Plot the actual data.
figure
plot(metersDistance, ... % X-axis
    meanRMS, ...
    'Marker', '*', 'MarkerSize', 4, ...
    'LineStyle', 'none','LineWidth',2, 'Color', [0 0 0]./255);
str_title = MatFile_str; ndx = strfind(str_title,'_'); str_title(ndx) = ' ';
title(str_title); set(gca,'FontSize',30); set(gca,'FontName','times')      
xlabel('Range (meters)');  ylabel('RMS Delay (ns)');
grid on;
ylim([fitvector(botgrab) fitvector(topgrab)])
%set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.

ax=gca;
ax.GridAlpha=0.3;

%jpegsavestring = ['RMSDelayPlotsScatter ' regexprep(MatFile_str,'_',' ','emptymatch')];
jpegsavestring = ['RMSDelayPlotsScatter'];
%jpegsavestring = savestring(1:end - 4);
jpegsavepath = [jpegsavepathbase '\' jpegsavestring];
jpegsavepath2 = jpegsavepath;

%%%Apparently a more exact way.
pause(0.0001)
hgexport(gcf, jpegsavepath2,  ...
    hgexport('factorystyle'), 'Format', 'jpeg');



%close all; 

%%%%%

% % % keyboard

%disp(['SAVE STRING: ', str_savefile]);
% save(str_savefileforus,'Strct_Metadata','t_rec','IQdata_Range_m','IQdata_Timing_cll', ...
%     'IQdata', 'pathlossfreqrec', 'pathlossavgrec','pathlossavgrec_dB','meanPathGain_dB','stdPathGain_dB',    '-v7.3');
%save(str_savefile,'Strct_Metadata','IQdata_CloudLocations_m_num','IQdata_Timing_cll','IQdata','-v7.3')

return
end

function [trueXX_interp,trueYY_interp,trueDistance,meanPathGain_dBp,meanPathGain_dBl] = distanceFromCenterOfCloud(Strct_Metadata,Txlocation,Rxlocation,Cloudlocations,meanPathGain)
% The locations are obained from the names of the TDMS files.
% locations should be a n by 2 array. The first column is the
% x coordinates, the second is the y coordinates.
    num = 101;    

    xPositions = Cloudlocations(:,1);
    trueXPosition = Rxlocation{1} - xPositions;
    trueXX = trueXPosition(1:sqrt(length(trueXPosition)):length(trueXPosition));
    trueXX_interp = linspace(min(trueXX),max(trueXX),num);
    
    yPositions = Cloudlocations(:,2);
    trueYPosition = Rxlocation{2} - yPositions;
    trueYY = trueYPosition(1:sqrt(length(trueYPosition)));
    trueYY_interp = linspace(min(trueYY),max(trueYY),num);
    
    [true_Xmesh,true_Ymesh] = meshgrid(trueXX,trueYY);
    [true_Xmesh_interp,true_Ymesh_interp] = meshgrid(trueXX_interp,trueYY_interp);
    
    trueDistance = sqrt((trueXPosition - Txlocation{1}).^2 + ...
                    (trueYPosition - Txlocation{2}).^2 + ...
                    (Rxlocation{3} - Txlocation{3}).^2);
  
    meanPathGain_reshape = reshape(meanPathGain,length(trueXX),length(trueYY));
    
    meanPathGain_interp_spline = interp2(true_Xmesh,true_Ymesh,meanPathGain_reshape, ...
                                        true_Xmesh_interp,true_Ymesh_interp,'spline');
    meanPathGain_dBp = 10.*log10(meanPathGain_interp_spline) - Strct_Metadata.ReceiverAntennaGain_dBi_num ...
                                 - Strct_Metadata.TransmitterAntennaGain_dBi_num;
    
    meanPathGain_interp_linear = interp2(true_Xmesh,true_Ymesh,meanPathGain_reshape, ...
                                         true_Xmesh_interp,true_Ymesh_interp,'linear');
    meanPathGain_dBl = 10.*log10(meanPathGain_interp_linear) - Strct_Metadata.ReceiverAntennaGain_dBi_num ...
                                 - Strct_Metadata.TransmitterAntennaGain_dBi_num;  
    
end


%%%RMS


function [trueXX_interp,trueYY_interp,trueDistance,rmsdata] = distanceFromCenterOfCloud_RMS(Txlocation,Rxlocation,Cloudlocations,rmsdata)
% The locations are obained from the names of the TDMS files.
% locations should be a n by 2 array. The first column is the
% x coordinates, the second is the y coordinates.
    num = 101;   
    
%     keyboard
    
    Txlocation = cell2mat(Txlocation);
    Rxlocation = cell2mat(Rxlocation);

    xPositions = Cloudlocations(:,1);
    trueXPosition = Rxlocation(1) - xPositions;
    trueXX = trueXPosition(1:sqrt(length(trueXPosition)):length(trueXPosition));
    trueXX_interp = linspace(min(trueXX),max(trueXX),num);
    
    yPositions = Cloudlocations(:,2);
    trueYPosition = Rxlocation(2) - yPositions;
    trueYY = trueYPosition(1:sqrt(length(trueYPosition)));
    trueYY_interp = linspace(min(trueYY),max(trueYY),num);
    
    [true_Xmesh,true_Ymesh] = meshgrid(trueXX,trueYY);
    [true_Xmesh_interp,true_Ymesh_interp] = meshgrid(trueXX_interp,trueYY_interp);
    
    trueDistance = sqrt((trueXPosition - Txlocation(1)).^2 + ...
                    (trueYPosition - Txlocation(2)).^2 + ...
                    (Rxlocation(3) - Txlocation(3)).^2);
  
    meanPathGain_reshape = reshape(rmsdata,length(trueXX),length(trueYY));
    
    rmsdata = interp2(true_Xmesh,true_Ymesh,meanPathGain_reshape, ...
                                        true_Xmesh_interp,true_Ymesh_interp,'spline');
    rmsdata = (rmsdata) ;
    
% % %     meanPathGain_interp_linear = interp2(true_Xmesh,true_Ymesh,meanPathGain_reshape, ...
% % %                                          true_Xmesh_interp,true_Ymesh_interp,'linear');
% % %     meanPathGain_dBl = 10.*log10(meanPathGain_interp_linear) - Strct_Metadata.ReceiverAntennaGain_dBi_num ...
% % %                                  - Strct_Metadata.TransmitterAntennaGain_dBi_num;  
    
end