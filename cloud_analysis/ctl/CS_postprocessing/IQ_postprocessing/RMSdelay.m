%Function that computes the average delay spread and RMS delay spread the
%imputs are PDP's  in column form and a timing vector in seconds the
%outputs are out1 = row vector of average delay time for each PDP
%out2  = row vector of Rms deleay for each PDP


function [timedelay,rmsdelay] = RMSdelay(PDP,t, pdp_threshhold)

t = t*1e9;% convert time to nanoseconds

%%%dynamic rescale yet?
dynamicrescale = 0;

% 

% taking pdp putting in log form not necassary help with visulization
Log_pdp = 10*log10(PDP);

%  "A" dB threshold

% pdp_threshhold = pdp_threshhold;
[~,pdpnumber] = size(PDP); % getting number of PDP'S put in to function

for indx = 1:pdpnumber
    
    threshold = max(Log_pdp(:,indx)) -pdp_threshhold; % threshold for each pdp %%%Max PDP - Thresh = Floor.
    
    %%%find the peak value of the given pdp
    peakdB = max(Log_pdp(:, indx) );
    %%%find the noise floor - 
    %%%%%noisefloordBapprox = 10*log10 (  mean(10.^(Log_pdp(8000:8100, indx)./10) ) ) ; %%%go through a chunk of the pdp VERY late in time. 
    
    noisefloordBapprox = max( Log_pdp(1500:1600, indx) ) ; %much more conservative estimate, by finding the MAX in this region. 
    
                    %%%this will approximate the long-term noise floor. 
    %%%Now, if the dB down (pdp_threshold) is less than this difference,
    %%%reduce the threshhold
    deltadb = peakdB - noisefloordBapprox;
    
    
    
     if deltadb < pdp_threshhold
         
         %%%dynamically rescale the pdp threshhold....
         tempthreshhold = deltadb - 3 ; %find the ACTUAL difference, and add 3dB for safety...
         
         threshold = max(Log_pdp(:,indx)) - tempthreshhold ; 
         
         dynamicrescale = 1;
         
         
         
         
     end
         
    
        
    
    %%%above threshhold
    indices = find(Log_pdp(:,indx)> threshold) ;
    endset = find(indices > 1600);
    indices(endset) = [];
    
    temp{indx} = indices ; % find indexes above threshold
    
    
    %     
    
end

% computing Average time delay equation 2b from ITU-R

for indx = 1:pdpnumber
    
%     
    %finding first peak (does both LOS and non-LOS)
    temp2 = Log_pdp(temp{indx}, indx );
    
    if size(temp2, 1) == 0
        continue;
    end
    
    [~,locs] = findpeaks(temp2);
    var = temp{indx};
    Tm = t(var(locs(1))); % time of first peak
    
    
    T = t(temp{indx});% times of points above threshold
    
    lpwr = PDP(temp{indx},indx); % power values above threshold
    
    Td(indx) = sum(T.*lpwr)./sum(lpwr)-Tm; % average time delay
    
    
    % computing RMS delay spread equation 4b from ITU-R
    
    Tm_vector = ones(size(T))*Tm; % so we dont have to use extra for loop
    Td_vector = ones(size(T))*Td(indx);
    
    S(indx) = (sum((T-Tm_vector-Td_vector).^2.*lpwr)./(sum(lpwr)))^.5; % RMS Delay spread
    
end

timedelay = Td;
rmsdelay = S;

if dynamicrescale == 1
    
   
    fprintf('DYNAMICALLY RESCALED THE RMS THRESHHOLD \n')
    fprintf(['WAS : ' num2str(pdp_threshhold) ' dB \n']);
    fprintf(['NOW: ' num2str(tempthreshhold) ' dB \n'] );
    
    
end

% 

