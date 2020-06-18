% Solve for impulse response and pdps of a single record

function [irpdp, Rec_wfm]=analyzeRecord_refs(Rec_wfm, Tx_wfm, t_rec, norm_atten,freqvector,flag_FilterType,filter_IQdata)

irpdp = cell(1, 2);


freqdata = (fft(Rec_wfm)./(fft(Tx_wfm)));

cbbcell_freq = freqdata.*filter_IQdata;
cbbcell = ifft(freqdata.*filter_IQdata./sum(filter_IQdata)).*norm_atten.*length(filter_IQdata);

filteravg = sum(filter_IQdata)./length(filter_IQdata);
filteravgsq = sum(filter_IQdata.^2)./length(filter_IQdata);

df = abs(freqvector(2) - freqvector(1));
dtrec = sqrt(t_rec(2) - t_rec(1)); %assuming constant sampling rate..
bwscale = sqrt(abs(freqvector(end) - freqvector(1)));
alphascale = 1./sum(abs(filter_IQdata.*freqdata).^2);

pnormfactor = sqrt( filteravgsq/length(filter_IQdata)/df) / filteravg;

%Pdp = (real(cbbcell).^2+imag(cbbcell).^2);

pnormpdp = cbbcell./sqrt(pnormfactor);

%Applying scaling
if (flag_FilterType == 0)   %0 = no filter
    PdpScaling = cbbcell;
elseif (flag_FilterType ==1)   %1 = PNIdeal with Area Sum Scale
    PdpScaling = pnormpdp./sqrt(bwscale);
elseif (flag_FilterType == 2)   %2 = PNIdeal with Discrete Sum Scale
    PdpScaling = pnormpdp .* sqrt(dtrec);
elseif (flag_FilterType == 3)   %3 = PNIdeal with Max Scale
    error('Max Scaling is currently not implemented in AnalyzeRecord_vs3')
elseif (flag_FilterType == 4)   %3 = Alpha Scaling for Path Gain
    PdpScaling_freq = cbbcell_freq .* sqrt(alphascale);
    PdpScaling = ifft(PdpScaling_freq).*norm_atten.*length(filter_IQdata);
end

irpdp{:,1} = PdpScaling; 

irpdp{:,2} = 10*log10(abs(irpdp{:,1}));

total_power = 10.*log10(sum(abs(PdpScaling_freq).^2))

clear n;    
keyboard
end