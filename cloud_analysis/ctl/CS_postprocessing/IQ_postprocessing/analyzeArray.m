%% solve for impulse response and Pdp of all records in an array

function [rec_irpdp, Rx_wfms]= analyzeArray(ArrayRecs, numRecords, wordlength, Tx_wfm,t_arr,t_rec, norm_atten)
 
    rec_irpdp=cell(numRecords,6);
    Rx_wfms=zeros(wordlength,numRecords,'double');

    for n=1:numRecords
        
        rec_irpdp{n,1}=t_arr(n);
        rec_irpdp{n,2}=(t_rec);
        
        for m=1:wordlength
            Rx_wfms(m,n) = ArrayRecs(m,n,1) + 1i*ArrayRecs(m,n,2);
            
        end

        rec_irpdp{n,3} = ifft(fft(Rx_wfms(:,n))./(fft(Tx_wfm)))*norm_atten;
        rec_irpdp{n,4} = 20*log10(abs(rec_irpdp{n,3}));
        %rec_irpdp{n,5} = rec_irpdp{n,4}(rec_irpdp{n,4} < max(rec_irpdp{n,4})-2);
        shift=wordlength/2;
        rec_irpdp{n,6} = circshift(rec_irpdp{n,4},[shift 0]);
    end
    
clear n; clear m;    

end

