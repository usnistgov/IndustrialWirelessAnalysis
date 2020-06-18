% get records into something useful

function [Record]=getRecord(MeasuredData, groupN, recordN, wordlength, numRecords)
    Record=zeros(wordlength,2,'double');
    prefix=2*numRecords*(groupN-1)+groupN+2;
    for n=1:wordlength
        Record(n,1)=MeasuredData(prefix+2*(recordN-1)).Data(n);
        Record(n,2)=MeasuredData(prefix+2*(recordN-1)+1).Data(n);
    end
end
