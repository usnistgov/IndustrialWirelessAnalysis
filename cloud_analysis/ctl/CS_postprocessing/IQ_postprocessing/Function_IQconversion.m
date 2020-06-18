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
function Function_IQconversion(MeasurementPath,t_file, t_arr, t_rec, Collection, header, ...
    numFiles, Tx_wfm, normalize,flag_FilterType,filter_IQdata,Strct_Metadata, oatsflag, min_amplitude)

%1. Get useful parameters from Strct_Metadata
numArrays =  Strct_Metadata.NumberAcqusitions_num;          %Number of Acqusitions
num_acqnum = Strct_Metadata.NumberAcqusitionExcelInfo_num; num_acqnum = [num_acqnum;numArrays];
num_dist = Strct_Metadata.Range_m_num; num_dist = [num_dist num_dist(end)];


%%%Construct the save file information.
MatFile_str = Strct_Metadata.MatFile_str;                   %Save Filename
MatFile_titlestr = MatFile_str;                            %title plots
MatFile_titlestr = regexprep(MatFile_titlestr,'_',' ','emptymatch'); %Trim the Underscores, which tends to make plot titles look bad.

%%%Create the str_save file, which names the file that we save.
str_savefile = [MeasurementPath MatFile_str '_pp.mat'];

str_savefileforus = [MeasurementPath MatFile_str '_ppforNIST.mat'];
str_savefileforrick = [MeasurementPath MatFile_str '_pp.mat'];
str_savefileforus = str_savefileforrick;

%%%More formatting.
Strct_Metadata.IQdata_Timing_Header_str = ['Acqusition Index' 'Record Index' 'Acqusition Time(date and seconds)' 'Record Time (seconds)' 'Total Time (date and seconds)'];
Strct_Metadata.IQdata_Range_Header_str = ['Acqusition Index' 'Record Index' 'Range Distance(m)' 'X(m)' 'Y(m)' 'Z(m)'];
Strct_Metadata.IQdata_Header_str = ['IQ Data (real and complex'];

numR=length(t_arr);  wordL=length(t_rec);
formatOut=('dd-mmm-yyyy HH:MM:SS:FFF');
%%%


%Create Frequency Conversion
%%%This creates the set of frequency points that correspond with the
%%%channel sounder measured data, based on the sampling time (dt), and the
%%%number of time points that are used (this is gathered from
%%% size(filter_IQdata) - a longer time stream will INCREASE frequency
%%% resolution (decrease the frequency bin size df), 

if isfield(Strct_Metadata, 'Frequency_GHz_num')
    freqc = Strct_Metadata.Frequency_GHz_num;
    freqc = freqc*1e9;
    dt = t_rec(2) - t_rec(1); 
    df = 1/dt;
    freqvec = transpose(-df/2 : df/size(filter_IQdata, 1) : df/2 - 1) + freqc;
    freqvec = freqvec./1e9; %convert to units of GHz.
    freqvector = freqvec.*1e9;
    
    
end

%%%Pre-Indexing thesecells.
IQdata_Timing_cll = cell(Strct_Metadata.NumberAcqusitions_num*Strct_Metadata.NumberRecordperAcqusition_num, 5);
IQdata_Range_m(Strct_Metadata.NumberAcqusitions_num*Strct_Metadata.NumberRecordperAcqusition_num, 1:3) = 0;


%%%%Run over the array information, construct various timing data.
%%%Construct Position Data as well.

strt_a = 1; ttl_ndx_acq = []; %record_distance = [];
for i = 1:size(num_acqnum,1)
    for a=strt_a:numArrays
        r = 1:numR;
        index=numR*(a-1)+r;
        
        IQdata_Timing_cll(index,1) = {a};  
        IQdata_Timing_cll(index,2) = num2cell(r);
        IQdata_Timing_cll(index,3) = t_file(a,2);           %Number Representation of Time
        IQdata_Timing_cll(index,4) = num2cell(seconds(t_arr(r)));
        IQdata_Timing_cll(index,5) = num2cell(datestr(datenum(t_file{a,2})+seconds(t_arr(r)),'dd-mmm-yyyy HH:MM:SS.FFF'),2);
        
        IQdata_Range_m(index,1)=a.*ones(1,numR);  IQdata_Range_m(index,2)=r;              
        if (a <= num_acqnum(i))
            ttl_ndx_acq = [ttl_ndx_acq index];
        else
            strt_a = a;
            break;
        end
    end
    
    
    %%%If we have the logical flag == true, then at oats, and the 
    %%%Lin space is a discrete step function (i.e. we stand in place for
    %%%some amount of time)
    if oatsflag == true
        if i == 1
            num_dist2 = linspace(num_dist(1),num_dist(1),size(ttl_ndx_acq,2));
        else
            num_dist2 = linspace(num_dist(i-1),num_dist(i-1),size(ttl_ndx_acq,2));
        end
    end
    %%%If we are NOT at the oats, then we have a moving measurement, and we
    %%%can approximate the position data as a linear fit over the
    %%%checkpoints.
    
    if oatsflag == false
        if i == 1
            num_dist2 = linspace(num_dist(1),num_dist(1),size(ttl_ndx_acq,2));
        else
            num_dist2 = linspace(num_dist(i-1),num_dist(i),size(ttl_ndx_acq,2));
        end
    end
    
    
    IQdata_Range_m(ttl_ndx_acq,3) = num_dist2;              %Record Distance (m)
    clear num_dist2 ndx_acq ttl_ndx_acq index
    ttl_ndx_acq = [];
    
end

%%%%


%%%Construct a linear interpolation of the position. 
%%%This introduces a small discontinuity at every change in 'checkpoint'
%%%(row of the Rvec matrix) - this can be addressed through a linear
%%%interpolation. Cubic Spline would be another approach. 
mylininterp = zeros(numArrays*numR, 1);
r0 = cell2mat(Strct_Metadata.Tx_xyz_m_cll);
x0 = r0(1); y0 = r0(2); z0 = r0(3); %The TX position - assumed constant. 
rvec = [Strct_Metadata.Rx_xyz_m_cll{1}; Strct_Metadata.Rx_xyz_m_cll{2}  ; Strct_Metadata.Rx_xyz_m_cll{3} ]';
%%%Construct the reciever position matrix - this is R vector.

starter = num_acqnum(1); %%%We have some set of records where we are in a fixed position. 

starter = Strct_Metadata.NumberRecordperAcqusition_num * starter ;

%%%X, Y, Z coordinates stored as well. 
xrecord = mylininterp;
yrecord = mylininterp;
zrecord = mylininterp;
mylininterp(1:starter) = sum( (rvec(1,:) - r0).^2).^(1/2);

%Initialize starting values. 
xrecord(1:starter) = (-rvec(1,1) + x0).^1;
yrecord(1:starter) = (-rvec(1,2) + y0).^1;
zrecord(1:starter) = (-rvec(1,3) + z0).^1;
%%%Main Loop. 
for s = 1:size(rvec,1) - 1
    
    x1 = rvec(s, 1); x2 = rvec(s+1, 1);
    y1 = rvec(s, 2); y2 = rvec(s+1, 2);
    %unchanging z, but...
    z1 = rvec(s, 3); z2 = rvec(s+1, 3);
    
    
    %So, to do this, I want to do a linear interpolation of the x and y
    %distances....   
    if s == 1
        acq1number = num_acqnum(s) + 0;
    else
        acq1number = num_acqnum(s) ;
    end
    
    acq1number;
    acq2number = num_acqnum(s+1);
   
    
    from = (acq1number-0)*numR + 1;
    thru = acq2number*numR;
    
    linx = linspace(x1, x2, (thru - from ) + 1);
    liny = linspace(y1, y2, (thru - from ) + 1);
    linz = linspace(z1, z2, (thru - from ) + 1);
    
    lind = (  (linx - x0).^2 + (liny - y0).^2 + (linz - z0).^2 ).^(1/2);
    
    
    mylininterp (  (from:thru) ) = lind;
    
    xrecord(from:thru) = (-linx + x0).^1;
    yrecord(from:thru) = (-liny + y0).^1;
    zrecord(from:thru) = (-linz + z0).^1; 
    
    if oatsflag == 1
        xrecord(from:thru) = x1;
        yrecord(from:thru) = y1;
        zrecord(from:thru) = z1;
        mylininterp(from:thru) = sqrt(  (x1 - x0).^2 + (y1- y0).^2 + (z1 - z0).^2 );
    end
end

%%%Edge Case caught here.

if (thru + 1) < length(mylininterp)
    
    mylininterp(thru+1 : end) = sum( (rvec(end,:) - r0).^2).^(1/2);
    
    xrecord(thru+1 : end) = (-rvec(end,1) + x0).^1;
    yrecord(thru+1 : end) = (-rvec(end,2) + y0).^1;
    zrecord(thru+1 : end) = (-rvec(end,3) + z0).^1;
    
end
%%%%%

%%%The IQ data Range contianing the Euclid. Distance, and vector distances
IQdata_Range_m(:, 3) = mylininterp;
IQdata_Range_m(:, 4) = -xrecord;
IQdata_Range_m(:, 5) = -yrecord;
IQdata_Range_m(:, 6) = -zrecord;

thruvec = 1:length(xrecord); %%% A useful vector for a variety of plotting - just a linear index vector.


%%%Plot of the Range data.
figure; 
plot(thruvec, IQdata_Range_m(:, 3), 'k', thruvec, IQdata_Range_m(:, 4), 'r--', ...
    thruvec, IQdata_Range_m(:, 5), 'b--', thruvec, IQdata_Range_m(:, 6) );
legend('Euclid Distance', 'X', 'Y', 'Z');
xlabel('Record number'); ylabel('D');
grid on;
title({['Distance plot vs record'], [MatFile_titlestr]});
%%%


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


% keyboard;

for j=1:numFiles
 
    for c=1:Collection{j,4}{11,2}
        
        %indexc=(numR*Collection{1,4}{11,2}*(j-1))+(numR*(c-1));
        %%%%Original index method. Tedious to use, replaced with 'indexer'.
        
        for t=1:numR
            
            indexer = indexer + 1;
            
            dummyRec=getRecord(Collection{j,1},c,t,wordL,numR);
            
            %[dummyPdp,~]=analyzeRecord(dummyRec,Tx_wfm,t_rec, normalize,flag_FilterType,filter_IQdata);
            
            %%%Original Analyze_Record - Results, excepting 'pathlossarea',
            %%%will agree extremely well. 
            
% % % %             [dummyPdp,~, pathlossfreqorig, pathlossavgorig, pathlossareaorig]=analyzeRecord_vs3(dummyRec,Tx_wfm,t_rec, ...
% % % %                 normalize,freqvector,flag_FilterType,filter_IQdata);
% % % %             
% % % % 
% % % % 
% % % %             IQdata_orig(:, indexer) = dummyPdp{:, 1};
% % % %             pathlossfreqrecavgorig = pathlossfreqrecavgorig + pathlossfreqorig;
% % % %             pathlossavgrecorig(indexer) = pathlossavgorig;
% % % %             pathlossarearecorig(indexer) = pathlossareaorig;

            
            %%%Cleaned up AnalyzeRecord Formulation.            

            [IQ_data, CIR_Discretescaling_Bandlimit, FreqResponse_allfreqtmp, PathGain_avgfreqtmp, PathGain_avgtimetmp]= ...
            analyzeRecord_vs5(dummyRec,Tx_wfm,t_rec,normalize,freqvector,flag_FilterType,filter_IQdata, vectorgrab);        
        
            raw_IQdata(:,indexer)= IQ_data;
            IQdata_scaled(:, indexer) = CIR_Discretescaling_Bandlimit;
            FreqResponse_allfreq(:,indexer) = FreqResponse_allfreqtmp;            
            PathGain_avgfreq(indexer) = PathGain_avgfreqtmp;
            PathGain_avgtime(indexer) = PathGain_avgtimetmp;
            
            
            %%%%PLFreq is identical. 
            

            

        end
    end
    
end

% keyboard;

% % % % % % %%%Finalize some bits - 
% % % % % % %%%Average the PathlossFreqRecAvg Approach.
% % % % % % pathlossfreqrecavg = pathlossfreqrecavg/(Strct_Metadata.NumberAcqusitions_num*Strct_Metadata.NumberRecordperAcqusition_num);
% % % % % % %%%FFTShift to be a little more readable. 
% % % % % % pathlossfreqrecavg = fftshift(pathlossfreqrecavg);
% % % % % % %%%%


%%%

figure; plotyy(thruvec, IQdata_Range_m(:, 3),  thruvec,10*log10(PathGain_avgfreq)); legend('Distance', 'Avg PL');
xlabel('Record'); ylabel('PLavg (dB)'); title({['ORIGINAL ACQ METHOD'],  [MatFile_titlestr]});

%%%%Plot various PDP from different records.
plot_color = {[0 51 51]./255, [120 100 60]./255, [120 100 20]./255};
ttl_figure = strtrim(MatFile_str);
ndx_figure = strfind(ttl_figure,'_'); ttl_figure(ndx_figure) = ' ';
figure; plot(IQdata_Range_m(:,3)); xlabel('Index of Record'); ylabel('Range (m)');
title(ttl_figure);

figure;
i1 = 100; AA1 = raw_IQdata(:,i1);
i2 = round(size(raw_IQdata,2)/2); AA2 = raw_IQdata(:,i2);
i3 = round(size(raw_IQdata,2)-100); AA3 = raw_IQdata(:,i3);
i4 = round(size(raw_IQdata,2) - 500); AA4 = raw_IQdata(:, i4);

plot(t_rec.*1e9,10.*log10(((real(AA1).^2+imag(AA1).^2))),'color',[0 0.85 0],'linewidth',2); hold on;
plot(t_rec.*1e9,10.*log10(((real(AA2).^2+imag(AA2).^2))),'color',[0 0.5 0],'linewidth',2);
plot(t_rec.*1e9,10.*log10(((real(AA3).^2+imag(AA3).^2))),'color',[0 0.15 0],'linewidth',2);
plot(t_rec.*1e9, ( 10.*log10(((real(AA4).^2+imag(AA4).^2))) ),'color',[0 0.15 0],'linewidth',2);
title(ttl_figure); grid on;
xlabel('Time (ns)'); ylabel('PDP (dB)'); set(gca,'FontSize',20); set(gca,'FontName','arial')
legend('Record = 100', ['Record = ' num2str(i2)], ['Record = ' num2str(i3)], ['Record = ' num2str(i4)]);
%axis([0 500 -180 -60])

%%%%Plot Path loss
thruvec = 1:length(PathGain_avgtime);
figure; hold on;
plot(IQdata_Range_m(:, 3), PathGain_avgtime,'o','color',[0 0.15 0],'linewidth',2)
set(gca,'FontSize',20); set(gca,'FontName','arial')
title(ttl_figure); ylabel('Path Gain (dB)'); xlabel('Range (m)')

figure;
plot(t_rec.*1e9,10.*log10(((real(AA1).^2+imag(AA1).^2))),'color',[0 0.85 0],'linewidth',2); hold on;
plot(t_rec.*1e9,10.*log10(((real(AA2).^2+imag(AA2).^2))),'color',[0 0.5 0],'linewidth',2);
plot(t_rec.*1e9,10.*log10(((real(AA3).^2+imag(AA3).^2))),'color',[0 0.15 0],'linewidth',2);
plot(t_rec.*1e9, ( 10.*log10(((real(AA4).^2+imag(AA4).^2))) ),'color',[0 0.15 0],'linewidth',2);

title(ttl_figure); grid on;
xlabel('Time (ns)'); ylabel('PDP (dB)'); set(gca,'FontSize',12); set(gca,'FontName','arial')
legend('Record = 100', ['Record = ' num2str(i2)], ['Record = ' num2str(i3)], ['Record = ' num2str(i4)]);
%axis([0 500 -180 -60])




[~, rmsvalues] = RMSdelay( abs(IQdata_scaled).^2 ,t_rec, 30 ) ;  %threshhold in db - timevector T

figure; plot(rmsvalues);
title('RMS Values');
ylabel('RMS Spread (ns)');
grid on;



%%%Path Setter.

MeasurementPath = 'Q:\NIST_Projects\SmartManufacturingData\Measurement Campaign Data\Processed Data\matfiles\';
str_savefileforrick = [MeasurementPath MatFile_str '_pp.mat'];

%%%This was a modified save file that was specific to a certain file being
%%%saved for Rick. 

% save(str_savefileforrick,'Strct_Metadata','IQdata_Range_m','IQdata_Timing_cll', ...
%      'IQdata','-v7.3');

 keyboard;
 
%%%Contains the full set of data.  

            save(str_savefileforrick,'Strct_Metadata','t_rec','IQdata_Range_m','IQdata_Timing_cll', ...
   'IQdata_scaled','FreqResponse_allfreq','PathGain_avgfreq', 'PathGain_avgtime', 'mylininterp',    '-v7.3');



 %save(str_savefileforus,'Strct_Metadata','t_rec','IQdata_Range_m','IQdata_Timing_cll', ...
 %    'IQdata', 'pathlossfreqrec', 'pathlossavgrec', 'mylininterp',    '-v7.3');
 
return





