function [header]=getHeader(IQ_File,groupNames,chanNames,refYN)
clear header;

header=cell(14,2);
for z=1:3
    header{z,1}=IQ_File.Data.Root.Property(z).Name;
    header{z,2}=IQ_File.Data.Root.Property(z).Value{1};
end
header(4,:)={'Date/Time',IQ_File.Data.Root.Property(4).Value};
for z=5:10
    header{z,1}=IQ_File.Data.Root.Property(z).Name;
    header{z,2}=IQ_File.Data.Root.Property(z).Value;
end
header(11,:)={'Number of Arrays',length(groupNames{1})};
header(12,:)={'Number of Records',(length(chanNames{1})/length(groupNames{1}))/2};
header(13,:)={'word length',IQ_File.Data.MeasuredData(3).Total_Samples};
header(14,:)={'WF increment',IQ_File.Data.MeasuredData(3).Property(1,3).Value};
if refYN==1
    header(15,:)={'Acquisition Delay (sec)',IQ_File.Data.MeasuredData(2).Property(3).Value};
    header(16,:)={'Max VST Input Level (dBm)',IQ_File.Data.MeasuredData(2).Property(4).Value};
end

return;

