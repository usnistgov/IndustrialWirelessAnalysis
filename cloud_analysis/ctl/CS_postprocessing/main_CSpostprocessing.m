%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%This main program provides a structure for the postporcessors of channel
%sounders.  These
% Function: main_CSpostprocessing
% Author: Jeanne Quimby - 03/17/2016
%         Alexandra Curtin
%         David Novotny
%--------------------------------------------------------------------------
%Inputs and Outputs for main_CSpostprocessing
%Inputs:
%
%Outputs
%   * num_distance
%--------------------------------------------------------------------------
%Function calls
%    none
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
close all; clear variables;

set(groot,'defaultAxesColorOrder',[31 120 180; 167 207 227; 51 160 44; 178 223 138; 227 26 28; 252 177 176; 255 127 0; 253 191 111]/255,'defaultAxesLineStyleOrder','-|--|:')
%0.0 Add matlab paths required for analysis
addpath('METADATA_postprocessing','IQ_postprocessing','PDP_postprocessing')

%0.1. Flags
%0.1.1 Perform both IQ and PDP postprocessing
flag_IQPDP = 1;       %0 = perform only IQ, 1 = perform both IQ and PDP, 2 = perform only PDP, 3 = check spreadsheet files are correct, 4 = look at references

min_amplitude = 0.1; %This limits the frequencies that are used to only be those that have |PNfilt|>Min Amp. Minamp = 0 uses all data. 0.1 is probably good.


%All filters currently perform area sum scaling on PN ideal scale.
flag_FilterType = 4;  %0 = no filter, 1 = PNIdeal with Area Sum Scale, 2 = PNIdeal with Discrete Sum Scale, 3 = PNIdeal with Max Scale, 4 = Area Sum Scaling
flag_CACSM = 1;        %0 = CAC measurements analysis = Pathloss, 1 = Smartmanufacturing Analysis  = Scaled IQ
cll_attentype{1} = 0; %0 - Assume a constant attenuation value. 1 - Assume a frequency dependent attenuation
if (cll_attentype{1} == 1)
    cll_attentype{2} = 'Q:\public\CAC_NIST\Analysis\VNA\Postprocessing\20160915\Analysis\B2B_50db\b2b50db.txt'; %50DB Atten .txt (IN DB FORM)
    cll_attentype{3} = 'Q:\public\CAC_NIST\Analysis\VNA\Postprocessing\20160915\Analysis\B2B_30dB\b2b30db.txt'; %30db Atten .txt (DB FORM)
else
    cll_attentype{2} = ''; cll_attentype{3} = '';
end

%that will be used to get better results. %NOTE - THIS NEEDS TO BE KEPT
%CONSISTENT WITH THE ATTENUATION USED IN THE B2B MEASUREMENT, OTHERWISE IT
%IWLL BE ERRONEOUS. SHOULD DO IN THE SPREADSHEET SOMEHOW

%str_fileexcelwholename  = 'Q:\public\CAC_NIST\Analysis\ChannelSounder\processeddata\CAC_MetaData.xlsx';

if (flag_CACSM == 0)
    
    str_fileexcelwholename = 'Q:\public\CAC_NIST\Analysis\Lab_1207_Measurements\ChannelSounder\processeddata\CAC_MetaData__vs2_copy.xlsx';
elseif (flag_CACSM == 1)
    %str_fileexcelwholename = 'Q:\NIST_Projects\SmartManufacturingData\Measurement Campaign Data\Processed Data\SmartManufacturing_MetaData2.xlsx';
    
    %str_fileexcelwholename = 'Q:\NIST_Projects\SmartManufacturingData\Measurement Campaign Data\Processed Data\SmartManufacturing_MetaData_cloud.xlsx';
    str_fileexcelwholename = 'Q:\NIST_Projects\SmartManufacturingData\software\CS_postprocessing\SmartManufacturing_MetaData_cloud.xlsx';
    
else
    error('CAC or Smartmanufacturing analysis not checked');
end

str_filterpath = './CS_parameters\filter\filter_Pnideal.mat';
load(str_filterpath);  %This loads the filter_IQdata
filter_IQdata = filter_Pnideal; %Awkward coding as is - THIS SHOULD BE CHANGED IN THE FUTURE

%1.0 Provide locations of measurement files
if (flag_IQPDP == 0)
    
    [num_all_filenames,text_allfilenames] = xlsread(str_fileexcelwholename,'Title Sheet','F32:I203');
    for ndx_file = 1:size(text_allfilenames,1)
        
        IQ_pathway = text_allfilenames{ndx_file,3};
        if (strcmpi(IQ_pathway,'\') == 0); IQ_pathway = [IQ_pathway '\']; end
        str_IQdatafile = text_allfilenames{ndx_file,1};
        %%%DETERMINE IF OATS?
        oatsflag = regexp(str_IQdatafile, 'Oats');
        if sum(oatsflag) >= 1
            oatsflag = logical(true);
        else
            oatsflag = logical(false);
        end
        
        disp(['STARTING METADATA PROCESSING']);
        
        [Strct_Metadata] = main_METADATApostprocessing(str_fileexcelwholename,str_IQdatafile);
        
        
        disp(['STARTING IQ POSTPROCESSING']);
        
        main_IQpostprocessing(IQ_pathway,flag_FilterType,filter_IQdata,Strct_Metadata, oatsflag, cll_attentype, min_amplitude);
        
    end
    
elseif(flag_IQPDP == 1)
    
    [num_all_filenames,text_allfilenames] = xlsread(str_fileexcelwholename,'Title Sheet','F32:I300');
    for ndx_file = 1:size(text_allfilenames,1)
        %for ndx_file = 1:1
        
        IQ_pathway = text_allfilenames{ndx_file,3};
        if (strcmpi(IQ_pathway,'\') == 0); IQ_pathway = [IQ_pathway '\']; end
        str_IQdatafile = text_allfilenames{ndx_file,1};
        if (strcmp(str_IQdatafile,'') == 1)
            'Completion of data runs'
            break;
        end
        
        
        %%%DETERMINE IF OATS?
        oatsflag = regexp(str_IQdatafile, 'Oats');
        if sum(oatsflag) >= 1
            oatsflag = logical(true);
        else
            oatsflag = logical(false);
        end
        %%%
        
        
        disp(['STARTING METADATA PROCESSING']);
        if (flag_CACSM == 0)
            [Strct_Metadata] = main_METADATApostprocessing_CAC(str_fileexcelwholename,str_IQdatafile);
        elseif (flag_CACSM == 1)
            [Strct_Metadata] = main_METADATApostprocessing(str_fileexcelwholename,str_IQdatafile);
        end
        disp(['STARTING IQ POSTPROCESSING']);
        
        main_IQpostprocessing(IQ_pathway,flag_FilterType,filter_Pnideal,Strct_Metadata,oatsflag, cll_attentype, min_amplitude);
    end
elseif(flag_IQPDP == 2)
    
    [num_all_filenames,text_allfilenames] = xlsread(str_fileexcelwholename,'Title Sheet','F32:I300');
    for ndx_file = 1:size(text_allfilenames,1)
        
        
        IQ_pathway = text_allfilenames{ndx_file,3};
        if (strcmpi(IQ_pathway,'\') == 0); IQ_pathway = [IQ_pathway '\']; end
        str_IQdatafile = text_allfilenames{ndx_file,1};
        
        str_IQdatamatfile = [IQ_pathway str_IQdatafile '_pp.mat'];
        flag_matfile = exist(str_IQdatamatfile);
        if (flag_matfile == 0)
            error([str_IQdatamatfile ' file does not exist'])
        end
        
        
        disp(['STARTING PDP POSTPROCESSING']);
        
        main_PDPpostprocessing_v4_quick(str_fileexcelwholename, str_IQdatamatfile, filter_Pnideal, [], [])
        
    end
elseif (flag_IQPDP == 3)  %check if files are correct.
    [num_all_filenames,text_allfilenames] = xlsread(str_fileexcelwholename,'Title Sheet','F32:I203');
    for ndx_file = 1:size(text_allfilenames,1)
        num_all_filenames(ndx_file)
        B2Bfile = text_allfilenames{ndx_file,2};
        IQ_pathway = text_allfilenames{ndx_file,3};
        if (strcmpi(IQ_pathway,'\') == 0); IQ_pathway = [IQ_pathway '\']; end
        if (strcmpi(B2Bfile(end-4:end),'.tdms') == 0); B2Bfile = [B2Bfile '.tdms']; end
        existDirectory = dir([IQ_pathway B2Bfile]);
        if (isempty(existDirectory))
            error([IQ_pathway B2Bfile ' does not exist']);
        end
    end
elseif (flag_IQPDP == 4)  %Look at multiple references
    
    [num_all_filenames,text_allfilenames] = xlsread(str_fileexcelwholename,'Title Sheet','F32:I203');
    [file_allinfonum,file_allinfotext] = xlsread(str_fileexcelwholename,'Test Runs Information','A1:T205');
    for ndx_file = 1:size(text_allfilenames,1)
        %1 Get pathways
        IQ_pathway = text_allfilenames{ndx_file,3};
        if (strcmpi(IQ_pathway,'\') == 0); IQ_pathway = [IQ_pathway '\']; end
        str_IQdatafile = text_allfilenames{ndx_file,1};
        if (strcmp(str_IQdatafile,'') == 1)
            'Completion of data runs'
            break;
        end
        
        Function_plotreferencefiles(IQ_pathway,filter_Pnideal);
        clear ndx_strct IQ_pathway str_IQdatafile
        
    end
end

