%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%This main program post processes .tdms files from labview channel sounder
%measurements.  Further documentation of the software can be found in Kate
%Remeley's report.
%Assumptions:
%   (1) multiple files are being processed
%   (2) the reference file is located in the same folder as the measurement
%   files
% Function: main_IQpostprocessing
% Author: Jeanne Quimby - 3/17/2016
%         David Novotny
%         Alexandra Curtin
%--------------------------------------------------------------------------
%Inputs and Outputs for main_ChannelSounderPostProcessing
%Inputs:
%   * none
%Outputs
%   * text file
%--------------------------------------------------------------------------
%Function calls
%    1) convertTDMS - converts labview TDMS files to matlab cells
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [] = main_IQpostprocessing(IQ_pathway,flag_FilterType,filter_IQdata,Strct_Metadata, oatsflag, cll_attentype, min_amplitude)

% 1.0 Reference File and Attenuation conversion from TDMS format and validation
%1.1 import data, read in headers and other useful numbers
ReferenceFile = [IQ_pathway Strct_Metadata.ReferenceFile_str];
if (strcmpi(ReferenceFile(end-4:end),'.tdms') == 0); ReferenceFile = [ReferenceFile '.tdms']; end
MeasurementPath = IQ_pathway;

%2 Conversion of Measurement Data
%2.0 Attenuation
atten = Strct_Metadata.Attenuation_dB_num;
if (cll_attentype{1} == 0)
    norm_atten=10^(-abs(atten)/20);
    cll_attentype{2} = atten;
    cll_attentype{3} = norm_atten;
elseif (cll_attentype{1} == 1)
    %
    % % %     %Figure out somehow how to load the vector attenuation.
    if (atten == 50)
        b2b_vec_atten = load(cll_attentype{2});
    elseif (atten == 30)
        b2b_vec_atten = load(cll_attentype{3});
    else
        error('Attenuation type is either 50 or 30, line 49 main_IQpostprocessing');
    end
    
    % % % % %     VNA and NIST CS are going to use very, very slightly different
    % % % % %     frequency vectors (probably). We can choose to assume they are close
    % % % % %     enough (this is likely fine), or do a cubic spline interpolation to
    % % % % %     attempt to sample A(w) at the correct frequency values. So
    
    %recreate frequency vector.
    %Create Frequency Conversion
    
    if isfield(Strct_Metadata, 'Frequency_GHz_num')
        freqc = Strct_Metadata.Frequency_GHz_num;
        freqc = freqc*1e9;
        %dt = t_rec(2) - t_rec(1);
        dt = 1/(Strct_Metadata.SampleRate_MHz_num*1e6);
        df = 1/dt;
        freqvec = transpose(-df/2 : df/size(filter_IQdata, 1) : df/2 - 1) + freqc;
        freqvec = freqvec./1e9; %convert to units of GHz.
        freqvector = freqvec.*1e9;
    else
        error('No Frequency Data for the frequency vector constructor, line 57, main_IQ_postprocessing');
        
    end
    
    %%%%%Find starting index....
    [~, a2l] = min (  abs(b2b_vec_atten(:, 1) - freqvec(1) ) ) ;
    %So we start at THIS value.
    %Go to where?
    [~, a2u] = min ( abs (b2b_vec_atten(:, 1) - freqvec(end) ) ) ;
    %Over whatever band we are doing.
    
    vna_atten = b2b_vec_atten(a2l:a2u, 2);
    vna_freqs = b2b_vec_atten(a2l:a2u + 0, 1);
    if length(vna_atten) < length(freqvec) %awkward fix...?
        vna_atten = b2b_vec_atten(a2l:a2u + 1, 2);
        vna_freqs = b2b_vec_atten(a2l:a2u + 1, 1);
        
    end
    %Convert to linear -
    vna_atten = 10.^( - abs(vna_atten)./20 );
    
    %CUBIC INTERPOLATION POSSIBILITY
    vnainterp = spline(vna_freqs, vna_atten, freqvec) ;
    %vnainterp = interp1(vna_freqs, vna_atten, freqvec, 'linear');
    %So now assign the b2b atten to this...
    norm_atten = vnainterp;
else
    error('Attenuation Type is incorrect, Line 40 in main_IQpostprocessing');
end
% keyboard

%2.1 Convert Reference data and header to matlab format
[Tx_File,~,Tx_chan,Tx_group] = convertTDMS(0,ReferenceFile);  %Converting Reference File to matlab cell

if size(regexp(ReferenceFile, 'cloud'), 1) == 0
    
   fprintf(['Converting: ' Strct_Metadata.ReferenceFile_str '\n']);
   
end

refmeta=getHeader(Tx_File,Tx_group,Tx_chan,1);
Ref_rec=getRecord(Tx_File.Data.MeasuredData,1,1,refmeta{13,2},refmeta{12,2});

%2.2 Average the Reference data
[Ref_allrecords, avgRef] =avgRecords(Tx_File.Data.MeasuredData, 1, refmeta{12,2}, refmeta{13,2});
avgRefWfm = avgRef(:,1)+ 1j*avgRef(:,2);
norm_avgRefWfm = avgRefWfm; %norm_avgRefWfm = avgRefWfm./sum(abs(avgRefWfm));
Ref_record1 =  Ref_allrecords(:,1,1) + 1j.*Ref_allrecords(:,1,2);
Ref_record2 =  Ref_allrecords(:,2,1) + 1j.*Ref_allrecords(:,2,2);
Ref_recordend =  Ref_allrecords(:,end,1) + 1j.*Ref_allrecords(:,end,2);

%2.3 Get reference timing only
[~,~,reft_rec]=createTiming(Tx_File.Data.MeasuredData,refmeta,1);

%2.8 Plot pdp from reference file

% 3.0 Post Processing of Received Data using Reference File
%--------------------------------------------------------------------------
% A header is: (1)name, (2)MeasNotes, (3)PNcode File, (4)Date/Time,
%(5)Carrier Freq, (6)Base PN code length, (7)PN Oversample, (8)Codeword length,
%(9)Acquire every N codewords, (10)Acquisitions (acquisitions) per file,
%(11)Number of acquisitions (pulls from group length this time, should match
%(12)Number of records, (13)Wordlength (this is Base word length * oversampling rate),
%(14)WF_increment (time between IQ data points), (15)Acquisition Delay between acquisitions,
%(16) Max VST Input Level.
%--------------------------------------------------------------------------
%3.1 Acquire Received Data File to convert
%[RecFiles,RecPaths]=uigetfile('*.tdms','Please choose all the *.tdms files that you need from this collection: ','MultiSelect','on',MeasurementPath);
MeasurementDirectory = dir([MeasurementPath '*.tdms']);
cll_MeasurementDirectory = struct2cell(MeasurementDirectory);
ndx_strct = cellfun(@(s) isempty(strfind(ReferenceFile, s)), cll_MeasurementDirectory(1,:));
ndx_meas = find(ndx_strct == 1);
IQ_datafile = sort(cll_MeasurementDirectory(1,ndx_meas));
numArrays=0; numFiles = size(IQ_datafile,2);

%%%Run through the set of detected TDMS file. 
for i=1:numFiles

    IQ_datafile{i};
    [IQ_dummy,~,Collection{i,3},Collection{i,2}]=convertTDMS(0,[MeasurementPath IQ_datafile{i}]);
    Collection{i,1}=IQ_dummy.Data.MeasuredData;
    Collection{i,4}=getHeader(IQ_dummy,Collection{i,2},Collection{i,3},1);
    numArrays=numArrays+Collection{i,4}{11,2};
    
end
%%%

%%%Script to create associated timing information for all the files.
[t_files, t_arr, t_rec]=createTiming(Collection, Collection{1,4},numFiles);
%%%

Y=['You have ',num2str(numFiles),' files, containing a total of ',num2str(numArrays),' acquisitions and ',num2str(Collection{1,4}{12,2}),' records per acquisition'];
disp(Y);


%4.0  Save Files and post processor
Strct_Metadata.NumberAcqusitions_num = numArrays;
Strct_Metadata.NumberFiles_num = numFiles;
Strct_Metadata.NumberRecordperAcqusition_num = Collection{1,4}{12,2};
%%%Cloud String check.
if (isempty(strfind(Strct_Metadata.MatFile_str,'Cloud')) == 0)
    [IQdata_CloudLocations_x, IQdata_CloudLocations_y] = getCloudLocations( IQ_pathway );
    IQdata_CloudLocations = [IQdata_CloudLocations_x IQdata_CloudLocations_y];
    Function_IQconversion_cloud(MeasurementPath,t_files, t_arr, t_rec, Collection, Collection{1,4}, ...
        numFiles, norm_avgRefWfm, norm_atten,filter_IQdata, Strct_Metadata,IQdata_CloudLocations, min_amplitude);
else
    %%%if NOT a cloud measurement, run through this processor. 
    Function_IQconversion(MeasurementPath,t_files, t_arr, t_rec, Collection, Collection{1,4}, ...
        numFiles, norm_avgRefWfm, norm_atten,flag_FilterType,filter_IQdata,Strct_Metadata, oatsflag, min_amplitude);
end
%%%


return

end

function [new_xCoordinates, new_yCoordinates] = getCloudLocations(folderpath)

files = dir([folderpath '/*.tdms']);

coordinatesArray = [];

for index = 1:size(files,1)
    filename = files(index).name;
    coordinatesArray = [coordinatesArray; fetchCoordinates(filename)];
end

coordinatesArray = coordinatesArray * (1/4000) * (1/39.3701);

if size(coordinatesArray, 2) == 2
    xCoordinates = coordinatesArray(:,1);
    yCoordinates = coordinatesArray(:,2);
else
    xCoordinates = [];
    yCoordinates = [];
end

%Find Center Point
midx = sum(xCoordinates)./length(xCoordinates);
midy = sum(yCoordinates)./length(yCoordinates);

new_xCoordinates = xCoordinates - midx;
new_yCoordinates = yCoordinates - midy;

end

function coordinates = fetchCoordinates(filename)

coordinatesCell = regexp(filename, '[+-]\d{7}', 'match');

if length(coordinatesCell) ~= 2
    coordinates = [];
else
    xCoordinate = str2double(coordinatesCell{1});
    yCoordinate = str2double(coordinatesCell{2});
    coordinates = [xCoordinate, yCoordinate];
end

end

