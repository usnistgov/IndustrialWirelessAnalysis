% Solve for impulse response and pdps of a single record

%%%Updated on 4/27/2017 by JDiener - Fixed the Frequency Energy Difference
%%%issue.

function [IQ_data, CIR_Discretescaling_Bandlimit, FreqResponse_allfreq, PathGain_avgfreq, PathGain_avgtime]= ...
    analyzeRecord_vs5(recorded_data, Tx_b2b, t_rec, norm_atten, freqvector, flag_FilterType, filter_IQdata, vectorgrab)



%1.0 Read in data from measured data, record data and other useful items
IQ_data = recorded_data(:, 1) + 1i.*recorded_data(:, 2); %Construct the Complex IQ waveform.
df = abs(freqvector(2) - freqvector(1));
dtrec = sqrt(t_rec(2) - t_rec(1)); %assuming constant sampling rate..

%2.0 Convert measured data to frequency and divided by the b2b.
freqdata = (fft(IQ_data)./(fft(Tx_b2b))); 

%2.1 Remove DC shift
%2.1.1 feq data ONE is the non-shifted zero dc point - or, fftshift, puts it at 4095...
freqdata(1) = 1/2 * (freqdata(2) + freqdata(end) ) ; %DC shift averaging...

%3.0 Filtering
%Code Block dealing with the set of frequency domain limiting.
%This is accomplished through a logical mask - we concern ourselves ONLY
%with the frequenceis above some relative amplitude level.

logicalmask = zeros( length(filter_IQdata), 1);
logicalmask(vectorgrab) = 1; %The square filter is implemented as a logicalmask.
filter_IQdata_bandlimit = filter_IQdata.*logicalmask;

%4.0 The 'Simple' Scaling method - this is what the scaling in the MUF
%reduces to 
Nf = length(filter_IQdata);
filtsumsq =  sum(filter_IQdata_bandlimit.^2);
discretescaling = sqrt(Nf)/sqrt(filtsumsq);


%5.0 Complex Base Band (CBB)
cbbcell = ifft(freqdata.*filter_IQdata_bandlimit./sum(filter_IQdata_bandlimit).*norm_atten).*length(filter_IQdata_bandlimit);

%6.0 Computes PDPs
CIR_Discretescaling_Bandlimit = ifft(freqdata.*filter_IQdata_bandlimit.*norm_atten).*discretescaling;

%7.0 Computes Path Gain
%7.1 Path Gain from raw frequency domain data
%Note: Nf is due to scaling inherent to matlab fft

% % FreqResponse_allfreq = ( abs(freqdata.*filter_IQdata_bandlimit.*norm_atten) .^2).*discretescaling./Nf; 
%%%original, above, had Discrete scaling outside. Put inside. 

FreqResponse_allfreq = ( abs(freqdata.*filter_IQdata_bandlimit.*norm_atten .*discretescaling).^2 )./Nf;

%7.2 Path Gain from average frequency domain data
%%%This is essentially an estimate of the path loss, by averaging over the
%%%frequency domain energy. 'Alpha' is a scaling factor that, for a flat
%%%channel, undoes the power of the filter. 

PathGain_avgfreq = sum( abs(freqdata.*filter_IQdata_bandlimit.*norm_atten.*discretescaling).^2 )./Nf;

%7.3 Path Gain from time domain data
%%%The estimate of the Path loss of the channel, as obtained from taking
%%%the sum of the PDP. This uses the truncated PDP. 
PathGain_avgtime = sum(abs(CIR_Discretescaling_Bandlimit).^2);

% keyboard


 
end