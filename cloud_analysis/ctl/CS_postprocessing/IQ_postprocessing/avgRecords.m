%% Create an IQ and AmpPhase dataset that is an average of all the records in an array

function [pull, avgRec]=avgRecords(MeasuredData, arrayN, numRecords, wordlength)

pull=getArray(MeasuredData, arrayN, numRecords, wordlength);
avgRec=zeros(wordlength,4,'double');
 for b=1:wordlength
     for c=1:numRecords
         avgRec(b,1)=avgRec(b,1)+pull(b,c,1);
         avgRec(b,2)=avgRec(b,2)+pull(b,c,2);
         avgRec(b,3)=avgRec(b,3)+sqrt(pull(b,c,1)^2+pull(b,c,2)^2);
         %This means that data for which I=0 causes problems in
         %the phase average.  THis also suggests there are places where
         %I=Acos(phi)=0 but Q=Asin(phi)!=0 = A, so phi=pi/2
         if pull(b,c,1)~=0
             avgRec(b,4)=avgRec(b,4)+atan((pull(b,c,2))/(pull(b,c,1)));
         elseif pull(b,c,1)==0; avgRec(b,4)=avgRec(b,4)+pi/2;end
     end
     avgRec(b,1)=avgRec(b,1)/numRecords; %avg I
     avgRec(b,2)=avgRec(b,2)/numRecords; %avg Q
     avgRec(b,3)=avgRec(b,3)/numRecords; %avg Amp
     avgRec(b,4)=avgRec(b,4)/numRecords; %avg Phase
 end
 
end