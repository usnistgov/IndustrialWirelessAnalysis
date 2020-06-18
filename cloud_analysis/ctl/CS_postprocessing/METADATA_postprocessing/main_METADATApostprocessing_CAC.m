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
[file_allinfonum,file_allinfotext] = xlsread(str_fileexcelwholename,'Test Runs Information','A11:S205');

%1.1.1 Excel file location for measurement filename
ndx_strct = cellfun(@(s) isempty(strfind(str_IQdatafile, s)), file_allinfotext(:,2));
if (isempty(ndx_strct) == 1); error('Error in reading file information in Test Runs Information sheet.  Check there is a \ at the end of the pathname or if the file exists in Column A in Test Runs Information.'); end
ndx1 = find(ndx_strct == 0);

%2 Transmitter and Receiver Location Information
%2.1 Transmitter Location

tmp_numtx = file_allinfonum(ndx1,11);
if (tmp_numtx == 1 && ndx1 < 205)  %Transmitter 1 from Lab 1207
    [num_txdistanceinfo,text_txdistanceinfo] = xlsread(str_fileexcelwholename,'Test Runs Information','E8:G8');
elseif (tmp_numtx == 1 && ndx1 >= 205) 
    [num_txdistanceinfo,text_txdistanceinfo] = xlsread(str_fileexcelwholename,'Test Runs Information','E9:E9');
else    
    error('Transmitter Location is incorrect.  Goto line 38 in matlab code.')
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

return

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

