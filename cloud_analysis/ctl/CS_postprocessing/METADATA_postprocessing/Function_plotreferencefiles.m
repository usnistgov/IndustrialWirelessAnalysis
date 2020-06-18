%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%This main program plots multiple reference files.
%Assumptions:
%   (1) multiple files are being processed
%   (2) the reference file is located in the same folder as the measurement
%   files
% Function: main_plotreferencefiles
% Author: Jeanne Quimby - 3/17/2016
%         David Novotny
%         Alexandra Curtin
%--------------------------------------------------------------------------
%Inputs and Outputs for none
%Inputs:
%   * none
%Outputs
%   * none
%--------------------------------------------------------------------------
%Function calls
%    1) cnone
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [] = Function_plotreferencefiles(MeasurementPath,filter_IQdata)

%0 Get files for reference plots
ndxslash = strfind(MeasurementPath,'\'); ReferencePath = [MeasurementPath(1:ndxslash(end-1)) 'referencefiles\'];
ReferenceDirectory = dir([ReferencePath '*.tdms']);
cll_ReferenceDirectory = struct2cell(ReferenceDirectory);
 
%2
IQ_datafile = sort(cll_ReferenceDirectory(1,:));
numArrays=0; numFiles = size(IQ_datafile,2);
for i=1:numFiles
    IQ_datafile{i}
        
    %2.1 Convert Reference data and header to matlab format
    [Tx_File,~,Tx_chan,Tx_group] = convertTDMS(0,[ReferencePath IQ_datafile{i}]);  %Converting Reference File to matlab acell
    refmeta=getHeader(Tx_File,Tx_group,Tx_chan,1);
    Ref_rec=getRecord(Tx_File.Data.MeasuredData,1,1,refmeta{13,2},refmeta{12,2});
    
    %2.2 Average the Reference data
    [Ref_allrecords, avgRef] =avgRecords(Tx_File.Data.MeasuredData, 1, refmeta{12,2}, refmeta{13,2});
    Refave = avgRef(:,1)+j*avgRef(:,2); fftRefave = fft(Refave); fftRefave(1) = 0.5.*(fftRefave(2)+fftRefave(end));
        
    Ref_record1 =  Ref_allrecords(:,1,1) + j.*Ref_allrecords(:,1,2); fftRef1 = fft(Ref_record1); fftRef1(1) = 0.5.*(fftRef1(2)+fftRef1(end));
    Ref_record2 =  Ref_allrecords(:,2,1) + j.*Ref_allrecords(:,2,2); fftRef2 = fft(Ref_record2); fftRef2(1) = 0.5.*(fftRef2(2)+fftRef2(end));
    Ref_recordend =  Ref_allrecords(:,end,1) + j.*Ref_allrecords(:,end,2); fftRefend = fft(Ref_recordend); fftRefend(1) = 0.5.*(fftRefend(2)+fftRefend(end));
        
    Ref_peak1=ifft(fftRef1.*filter_IQdata./(fftRefave.*sum(filter_IQdata))).*length(filter_IQdata);
    Refpdp1 = 20*log10(abs(Ref_peak1));
    
    Ref_peak2=ifft(fftRef1.*filter_IQdata./(fftRefave.*sum(filter_IQdata))).*length(filter_IQdata);
    Refpdp2 = 20*log10(abs(Ref_peak2));
    Ref_peakend=ifft(fftRefend.*filter_IQdata./(fftRefave.*sum(filter_IQdata))).*length(filter_IQdata);
    Refpdpend = 20*log10(abs(Ref_peakend));
    
        
    %2.3 Get reference timing only
    [~,~,reft_rec]=createTiming(Tx_File.Data.MeasuredData,refmeta,1);
    
    %2.8 Plot pdp from reference file
    figure
    plot(reft_rec./1e-6,Refpdp1,'k',reft_rec./1e-6,Refpdp2,'b',reft_rec./1e-6,Refpdpend,'r');
    ttl_ReferenceFile = strtrim( IQ_datafile{i});
    ndx_ReferenceFile = strfind(ttl_ReferenceFile,'_'); ttl_ReferenceFile(ndx_ReferenceFile) = ' ';
    title(['Reference Peak for ' ttl_ReferenceFile])
    xlabel('Time (microseconds)'); ylabel('Magnitude (dB)')
    legend('Record 1 divided Ave Records','Record 2 divided Ave Records', ...
        'Last Record divided Ave Records')
    grid on;
    keyboard
    clear Ref_rec avgRef reft_rec Tx_File Tx_chan Tx_group
end


return