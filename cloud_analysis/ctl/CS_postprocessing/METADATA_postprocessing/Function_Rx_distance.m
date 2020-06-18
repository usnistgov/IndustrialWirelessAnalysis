%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%This function computes the Receiver distance based upon the spreadsheet true and
% measuremetn values of x and ys
% Function: Function_Rx_distance
% Author: Jeanne Quimby - 11/13/2015
%         David Novotny
%--------------------------------------------------------------------------
%Inputs and Outputs
%Inputs:
%   *
%Outputs
%   *
%--------------------------------------------------------------------------
%Function calls
%    none
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [num_acqnum,num_rx_x, num_rx_y,num_rx_z] = Function_Rx_distance(file_excelwholename,str_IQdatafile,distance_type)

%1.3 Receiver Location Information
[Rx_positionnum,Rx_positiontext] = xlsread(file_excelwholename,'Data_AcqLoc','A5:AV430'); % read data_acqloc excel sheet
[ndx1,ndx2] = find(    cellfun( 'length',   regexp(Rx_positiontext,str_IQdatafile)    ) == 1);
if (isempty(ndx1) == 1); error('Error in reading file information in Data_AcqLoc sheet.'); end
ndx1 = ndx1(1); ndx2 = ndx2(1);

%1.3.0 Determine Number of Acqusitions
num_acqnum(1) = Rx_positionnum(ndx1,ndx2+1);
tmp_num_acqnum = Rx_positionnum(:,ndx2+1);  
ndx_nan_acqnum = find(     isnan(tmp_num_acqnum) == 1); % finding where NAN is in the file
ndx2_nan_acqnum = find(ndx_nan_acqnum - ndx1 > 0); % at moment there are no nans beyond our acquisistions

if (isempty(ndx2_nan_acqnum) == 1)
    num_acqnum = Rx_positionnum(ndx1:end,ndx2+1); % number aquisistions
    num_turnnum = Rx_positionnum(ndx1:end,1); % checkpoints
    text_movement = Rx_positiontext(ndx1:end,ndx2+1); % type of movement
else
    num_acqnum = Rx_positionnum(ndx1:ndx_nan_acqnum(ndx2_nan_acqnum(1))-1,ndx2+1);
    num_turnnum = Rx_positionnum(ndx1:ndx_nan_acqnum(ndx2_nan_acqnum(1))-1,1);
    text_movement = Rx_positiontext(ndx1:ndx_nan_acqnum(ndx2_nan_acqnum(1))-1,ndx2+1);
end 

[Rx_checknum,Rx_checktext] = xlsread(file_excelwholename,distance_type,'A4:BI44');   %*****
%[Rx_checknum,Rx_checktext] =
%xlsread(file_excelwholename,distance_type,'A4:AQ44');  %with most cup runs

%1.3.2 Determine the Number of Acqusitions versus Number of Turns versus
num_rx_x = []; num_rx_y = []; num_rx_z = []; ii = 1;
for i = ndx1:size(num_acqnum,1)+ndx1-1 
    
%for i = ndx1:size(num_acqnum,1)+ndx1-2 %(MADE THIS CHANGE KEPT OVER INDEXING - Chad)

    if (strcmpi(Rx_positiontext(i,ndx2+1),'Stationary') == 1)
        ndx_check = find(Rx_checknum(:,6) == num_turnnum(i-ndx1+1)); ndx_check = ndx_check(1);
        ttl_num_acqnum(ii,1) = num_acqnum(i-ndx1+1);
        [num_rxdistanceinfo,text_rxdistanceinfo] = xlsread(file_excelwholename,distance_type,'C4:E38'); 
    elseif (strcmpi(Rx_positiontext(i,ndx2+1),'Ignore') == 1)
        continue
    elseif (strcmpi(Rx_positiontext(i,ndx2+1),'Walking') == 1)
        ndx_check = find(Rx_checknum(:,6) == num_turnnum(i-ndx1+1));
        if (size(ndx_check,1) == 1)
            ttl_num_acqnum(ii,1) = num_acqnum(i-ndx1+1);
        else
            keyboard
            ttl_num_acqnum(ii:ii+size(ndx_check,1)-1,1) = linspace(num_acqnum(i-ndx1+1),num_acqnum(i-ndx1+2)-1,size(ndx_check,1));
        end
        [num_rxdistanceinfo,text_rxdistanceinfo] = xlsread(file_excelwholename,distance_type,'C4:E38');
       
    elseif (strcmpi(Rx_positiontext(i,ndx2+1),'Inner Loop, Stationary') == 1)
        ndx_check = find(Rx_checknum(:,24) == num_turnnum(i-ndx1+1)); ndx_check = ndx_check(1);
        ttl_num_acqnum(ii,1) = num_acqnum(i-ndx1+1);
        [num_rxdistanceinfo,text_rxdistanceinfo] = xlsread(file_excelwholename,distance_type,'U4:W25');
        
    elseif (strcmpi(Rx_positiontext(i,ndx2+1),'Inner Loop, Walking') == 1)
        ndx_check = find(Rx_checknum(:,24) == num_turnnum(i-ndx1+1));
        if (size(ndx_check,1) == 1)
            ttl_num_acqnum(ii,1) = num_acqnum(i-ndx1+1);
        else
            ttl_num_acqnum(ii:ii+size(ndx_check,1)-1,1) = linspace(num_acqnum(i-ndx1+1),num_acqnum(i-ndx1+2)-1,size(ndx_check,1));
        end
        [num_rxdistanceinfo,text_rxdistanceinfo] = xlsread(file_excelwholename,distance_type,'U4:W25');

    elseif (strcmpi(Rx_positiontext(i,ndx2+1),'Outer Loop, Stationary') == 1)
        ndx_check = find(Rx_checknum(:,18) == num_turnnum(i-ndx1+1)); ndx_check = ndx_check(1);
        ttl_num_acqnum(ii,1) = num_acqnum(i-ndx1+1);
        [num_rxdistanceinfo,text_rxdistanceinfo] = xlsread(file_excelwholename,distance_type,'O4:Q29');
    elseif (strcmpi(Rx_positiontext(i,ndx2+1),'Outer Loop, Walking') == 1)
        ndx_check = find(Rx_checknum(:,18) == num_turnnum(i-ndx1+1));        
        if (size(ndx_check,1) == 1)
            ttl_num_acqnum(ii,1) = num_acqnum(i-ndx1+1);
        else
            ttl_num_acqnum(ii:ii+size(ndx_check,1)-1,1) = linspace(num_acqnum(i-ndx1+1),num_acqnum(i-ndx1+2)-1,size(ndx_check,1));
        end
        [num_rxdistanceinfo,text_rxdistanceinfo] = xlsread(file_excelwholename,distance_type,'O4:Q29');

    elseif (strcmpi(Rx_positiontext(i,ndx2+1),'CUP, Stationary') == 1)
        ndx_check = find(Rx_checknum(:,31) == num_turnnum(i-ndx1+1)); ndx_check = ndx_check(1);
        ttl_num_acqnum(ii,1) = num_acqnum(i-ndx1+1);
        [num_rxdistanceinfo,text_rxdistanceinfo] = xlsread(file_excelwholename,distance_type,'AG4:AI41');
    elseif (strcmpi(Rx_positiontext(i,ndx2+1),'CUP, Walking') == 1)
        ndx_check = find(Rx_checknum(:,31) == num_turnnum(i-ndx1+1));        
        if (size(ndx_check,1) == 1)
            ttl_num_acqnum(ii,1) = num_acqnum(i-ndx1+1);
        else
            ttl_num_acqnum(ii:ii+size(ndx_check,1)-1,1) = linspace(num_acqnum(i-ndx1+1),num_acqnum(i-ndx1+2)-1,size(ndx_check,1));
        end
        [num_rxdistanceinfo,text_rxdistanceinfo] = xlsread(file_excelwholename,distance_type,'AG4:AI41');
    elseif (strcmpi(Rx_positiontext(i,ndx2+1),'CUP OATS') == 1)
    
        ndx_check = find(Rx_checknum(:,37) == num_turnnum(i-ndx1+1));  
  
        if (size(ndx_check,1) == 1)
            ttl_num_acqnum(ii,1) = num_acqnum(i-ndx1+1);
        else
            ttl_num_acqnum(ii:ii+size(ndx_check,1)-1,1) = linspace(num_acqnum(i-ndx1+1),num_acqnum(i-ndx1+2)-1,size(ndx_check,1));
        end
        %[num_rxdistanceinfo,text_rxdistanceinfo] = xlsread(file_excelwholename,distance_type,'AM4:AO44');   %***** Correct for most
        [num_rxdistanceinfo,text_rxdistanceinfo] = xlsread(file_excelwholename,distance_type,'AM4:AO48');        
     elseif (strcmpi(Rx_positiontext(i,ndx2+1),'Lobby Walking') == 1)
        ndx_check = find(Rx_checknum(:,54) == num_turnnum(i-ndx1+1));
        if (size(ndx_check,1) == 1)
            ttl_num_acqnum(ii,1) = num_acqnum(i-ndx1+1);
        else
            ttl_num_acqnum(ii:ii+size(ndx_check,1)-1,1) = linspace(num_acqnum(i-ndx1+1),num_acqnum(i-ndx1+2)-1,size(ndx_check,1));
        end
        [num_rxdistanceinfo,text_rxdistanceinfo] = xlsread(file_excelwholename,distance_type,'AY4:BA7');

    else error(['Undefined Movement Type in Excel Spreadsheet, ' Rx_positiontext{i,ndx2+1}]);end

    %1.3.1 Determine the Receiver Locations
    
    num_rx_x = cat(1,num_rx_x,num_rxdistanceinfo(ndx_check,1));
    num_rx_y = cat(1,num_rx_y,num_rxdistanceinfo(ndx_check,2));
    num_rx_z = cat(1,num_rx_z,num_rxdistanceinfo(ndx_check,3));
   ii = ii+size(ndx_check,1);
   
end
num_acqnum = round(ttl_num_acqnum);

%keyboard;

return

