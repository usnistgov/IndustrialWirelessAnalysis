%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%This function computes the distance based upon the acqusition number in
%the .mat files
% Function: main_METADATA_postprocessing
% Author: Jeanne Quimby - 03/21/2016
%         David Novotny
%--------------------------------------------------------------------------
%Inputs and Outputs
%Inputs:
%   * num_acqusition file
%Outputs
%   * num_distance
%--------------------------------------------------------------------------
%Function calls
%    none
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [Strct_Metadata] = main_METADATApostprocessing(str_fileexcelwholename,str_IQdatafile)

%1.0 Compute Distance from excel spreadsheet and knowledge of transmitter and receiver location
%1.1 Load excel spreadsheet
[file_allinfonum,file_allinfotext] = xlsread(str_fileexcelwholename,'Test Runs Information','A13:S109'); % seperates numbers and strings in excel sheet

%1.1.1 Excel file location for measurement filename
ndx_strct = cellfun(@(s) isempty(strfind(str_IQdatafile, s)), file_allinfotext(:,2)); % looks in column 2 for  str_IQdatafile
if (isempty(ndx_strct) == 1); error('Error in reading file information in Test Runs Information sheet.  Check there is a \ at the end of the pathname or if the file exists in Column A in Test Runs Information.'); end
ndx1 = find(ndx_strct == 0);% finding index where str_IQdatafile is located

%2 Transmitter and Receiver Location Information
%2.1 Transmitter Location
% 
  

tmp_numtx = file_allinfonum(ndx1,11); 

% chad u have to update this according to your infor on the metadata excel
% sheet test runs info
% 'This could be causing problems'

 %CHANGING SOME < 19 -> <=  , and the >= 19 to > 19
if (tmp_numtx == 1 && ndx1 <= 19)  %Transmitter 1 from Gaithersburg Data
    [num_txdistanceinfo,text_txdistanceinfo] = xlsread(str_fileexcelwholename,'Test Runs Information','E8:G8');
elseif (tmp_numtx == 2 && ndx1 <= 19) %Transmitter 2 from Gaithersburg Data
    [num_txdistanceinfo,text_txdistanceinfo] = xlsread(str_fileexcelwholename,'Test Runs Information','E9:G9');
elseif (tmp_numtx == 3 && ndx1 <= 19) %TX 3 From Gaithesburg....
    [num_txdistanceinfo,text_txdistanceinfo] = xlsread(str_fileexcelwholename,'Test Runs Information','E10:G10');
elseif (tmp_numtx == 1 && ndx1 > 19 && ndx1 < 45) %Transmitter 1 from Automotative Assembly Data
    [num_txdistanceinfo,text_txdistanceinfo] = xlsread(str_fileexcelwholename,'Test Runs Information','E34:G34');
elseif (tmp_numtx == 2 && ndx1 > 19 && ndx1 < 45) %Transmitter 2 from Automotative Assembly Data
    [num_txdistanceinfo,text_txdistanceinfo] = xlsread(str_fileexcelwholename,'Test Runs Information','E35:G35');
elseif (tmp_numtx == 1 && ndx1 >= 45 &&ndx1 < 96) % Omnidirectional Transmitter at Oats
    [num_txdistanceinfo,text_txdistanceinfo] = xlsread(str_fileexcelwholename,'Test Runs Information','E55:G55');
elseif (tmp_numtx == 2 && ndx1 >= 45) %Transmitter 2 from CUP in Boulder, NIST
    [num_txdistanceinfo,text_txdistanceinfo] = xlsread(str_fileexcelwholename,'Test Runs Information','E56:G56');
elseif (tmp_numtx == 3 && ndx1 >= 45) %Transmitter 3 from CUP at Boulder, NIST
    [num_txdistanceinfo,text_txdistanceinfo] = xlsread(str_fileexcelwholename,'Test Runs Information','E57:G57');
elseif (tmp_numtx == 4 && ndx1 >= 45) % Horn Transmitter at Oats
    [num_txdistanceinfo,text_txdistanceinfo] = xlsread(str_fileexcelwholename,'Test Runs Information','E58:G58');
%% Lobby measurements with alex
elseif(tmp_numtx ==1 && ndx1 >= 96) %omni antenna in the lobby of nist
      [num_txdistanceinfo,text_txdistanceinfo] = xlsread(str_fileexcelwholename,'Test Runs Information','E103:G103')
   
    

else
    
    error('Transmitter Location is incorrect.  Goto line 38 in matlab code.')
end

%Quick NISTCS fix to Metadata processing...
if sum ( regexp(str_IQdatafile ,'NISTCS')) > 0 
   %Now look up data...
   if tmp_numtx == 1
   [num_txdistanceinfo,text_txdistanceinfo] = xlsread(str_fileexcelwholename,'Test Runs Information','E8:G8');
   elseif tmp_numtx == 2
      [num_txdistanceinfo,text_txdistanceinfo] = xlsread(str_fileexcelwholename,'Test Runs Information','E9:G9'); 
   end
       
end
num_tx_x = num_txdistanceinfo(1); num_tx_y = num_txdistanceinfo(2); num_tx_z = num_txdistanceinfo(3);

%2.2 Receiver Location
distance_type = 'Data_Physical_Location';

[num_acqnum,num_rx_x, num_rx_y, num_rx_z] = Function_Rx_distance(str_fileexcelwholename,str_IQdatafile,distance_type);

%2.3 Compute distance between transmitter and receiver
num_dist = sqrt((num_rx_x - num_tx_x).^2 + (num_rx_y - num_tx_y).^2 + (num_rx_z - num_tx_z).^2).';

%3. Return Cell Array of METADATA
Strct_Metadata.MatFile_str = file_allinfotext{ndx1,2};
Strct_Metadata.Frequency_GHz_num = file_allinfonum(ndx1,4);
Strct_Metadata.Location_str = [file_allinfotext{ndx1,3} ' at ' file_allinfotext{ndx1,5}];
Strct_Metadata.ReceiverAntenna_str = file_allinfotext{ndx1,6};
Strct_Metadata.ReceiverAntennaGain_dBi_num = file_allinfonum(ndx1,7);
Strct_Metadata.TransmitterAntenna_str = file_allinfotext{ndx1,8};
Strct_Metadata.TransmitterAntennaGain_dBi_num = file_allinfonum(ndx1,9);
Strct_Metadata.TransmitterPower_watts_num = file_allinfonum(ndx1,10);
Strct_Metadata.ReferenceFile_str = file_allinfotext{ndx1,13};
Strct_Metadata.Attenuation_dB_num = file_allinfonum(ndx1,14);
Strct_Metadata.PnChipRate_num = file_allinfonum(ndx1,15);
Strct_Metadata.BasePNCcodeLength_num = file_allinfonum(ndx1,16);
Strct_Metadata.PNOversample_num = file_allinfonum(ndx1,17);
Strct_Metadata.CodewordLength_num = file_allinfonum(ndx1,18);
Strct_Metadata.SampleRate_MHz_num = file_allinfonum(ndx1,19);
Strct_Metadata.Range_m_num = num_dist;
Strct_Metadata.Rx_xyz_m_cll = {num_rx_x.',num_rx_y.',num_rx_z.'};
Strct_Metadata.Tx_xyz_m_cll = {num_tx_x.',num_tx_y.',num_tx_z.'};
Strct_Metadata.NumberAcqusitionExcelInfo_num = num_acqnum;

 %keyboard;

return

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

