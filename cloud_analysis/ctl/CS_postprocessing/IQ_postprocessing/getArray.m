%%get an array of records out

function [myArray]=getArray(MeasuredData, groupN, numRecords,wordlength)
    clear prefix; clear Array;
    k=0;
    
    myArray=zeros(wordlength, numRecords,2,'double');
    prefix=2*numRecords*(groupN-1)+groupN+1;
    for n=1:(numRecords*2)
        if mod(n,2)==1
            k=k+1;
            for m=1:(wordlength)
                myArray(m,k,1)=MeasuredData(prefix+n).Data(m);
            end
        elseif mod(n,2)==0
            for p=1:(wordlength)
                myArray(p,k,2)=MeasuredData(prefix+n).Data(p);
            end
        end
    end
    clearvars n m p k;
end 