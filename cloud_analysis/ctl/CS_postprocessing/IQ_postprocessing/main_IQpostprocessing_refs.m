function [] = main_IQpostprocessing_refs(IQ_pathway,flag_FilterType,filter_IQdata,Strct_Metadata)

% 1.0 Reference File and Attenuation conversion from TDMS format and validation
%1.1 import data, read in headers and other useful numbers
ReferenceFile = [IQ_pathway Strct_Metadata.ReferenceFile_str];
if (strcmpi(ReferenceFile(end-4:end),'.tdms') == 0); ReferenceFile = [ReferenceFile '.tdms']; end
MeasurementPath = IQ_pathway;

%2 Conversion of Measurement Data
%2.0 Attenuation
atten = Strct_Metadata.Attenuation_dB_num;
norm_atten=10^(-abs(atten)/20);

%2.1 Convert Reference data and header to matlab format
[Tx_File,~,Tx_chan,Tx_group] = convertTDMS(0,ReferenceFile);  %Converting Reference File to matlab acell
refmeta=getHeader(Tx_File,Tx_group,Tx_chan,1);
Ref_rec=getRecord(Tx_File.Data.MeasuredData,1,1,refmeta{13,2},refmeta{12,2});

%2.2 Average the Reference data
[Ref_allrecords, avgRef] =avgRecords(Tx_File.Data.MeasuredData, 1, refmeta{12,2}, refmeta{13,2});
avgRefWfm = avgRef(:,1)+j*avgRef(:,2);
norm_avgRefWfm = avgRefWfm; %norm_avgRefWfm = avgRefWfm./sum(abs(avgRefWfm));
Ref_record1 =  Ref_allrecords(:,1,1) + j.*Ref_allrecords(:,1,2);
Ref_record2 =  Ref_allrecords(:,2,1) + j.*Ref_allrecords(:,2,2);
Ref_recordend =  Ref_allrecords(:,end,1) + j.*Ref_allrecords(:,end,2);

Ref_peak1=ifft(fft(Ref_record1)./(fft(avgRefWfm)));
Refpdp1 = 20*log10(abs(Ref_peak1));
Ref_peak2=ifft(fft(Ref_record2)./(fft(avgRefWfm)));
Refpdp2 = 20*log10(abs(Ref_peak2));
Ref_peakend=ifft(fft(Ref_recordend)./(fft(avgRefWfm)));
Refpdpend = 20*log10(abs(Ref_peakend));
Ref_peak21 = ifft(fft(Ref_record2)./(fft(Ref_record1)));
Refpdp21 = 20*log10(abs(Ref_peak21));
Ref_peakend1 = ifft(fft(Ref_recordend)./(fft(Ref_record1)));
Refpdpend1 = 20*log10(abs(Ref_peakend1));

%2.3 Get reference timing only
[~,~,reft_rec]=createTiming(Tx_File.Data.MeasuredData,refmeta,1);

%2.4 Plot pdp from reference file
figure
plot(reft_rec,Refpdp1,'k',reft_rec,Refpdp2,'b',reft_rec,Refpdpend,'r', ...
    reft_rec,Refpdp21,'c',reft_rec,Refpdpend1,'g');
ttl_ReferenceFile = strtrim(Strct_Metadata.ReferenceFile_str);
ndx_ReferenceFile = strfind(ttl_ReferenceFile,'_'); ttl_ReferenceFile(ndx_ReferenceFile) = ' ';
title(['Reference Peak for ' ttl_ReferenceFile])
xlabel('Time (s)'); ylabel('Reference File with Attenuation (dB)')
legend('Record 1 divided Ave Records','Record 2 divided Ave Records', ...
    'Last Record divided Ave Records','Record 2 divided Record 1','Last Record divided Record 1')
grid on;

%2.5 Determine Scaling Factor
%2.5.1 %Create Frequency Conversion
if isfield(Strct_Metadata, 'Frequency_GHz_num')
    freqc = Strct_Metadata.Frequency_GHz_num;
    freqc = freqc*1e9;
    dt = reft_rec(2) - reft_rec(1); df = 1/dt;
    freqvec = transpose(-df/2 : df/size(filter_IQdata, 1) : df/2 - 1) + freqc;
    freqvec = freqvec./1e9; %convert to units of GHz.
    freqvector = freqvec.*1e9;
    
end

%2.5.2 Convert Reference File
[irpdp, Rec_wfm]=analyzeRecord_refs(Ref_record1, norm_avgRefWfm, reft_rec, norm_atten,freqvector,flag_FilterType,filter_IQdata)



keyboard