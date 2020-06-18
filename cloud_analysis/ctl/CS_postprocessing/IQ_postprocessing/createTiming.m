%% Make some arrays about timing

function [t_file, t_arr,t_rec]=createTiming(MeasuredData, header, numFiles)
t_file=cell(header{11,2},2);

t_rec=zeros(header{13,2},1,'double');
for a=1:header{13,2}
    t_rec(a)=a*header{14,2};
end

wordtime=header{13,2}*header{14,2};  %I'm going to need to rewrite how this is defined
t_arr=zeros(header{12,2},1,'double');
for b=1:header{12,2}
    t_arr(b)=b*(wordtime*header{9,2}+header{15,2});
end

if numFiles==1
    d=1;
    if (iscell(MeasuredData) == 1)
        for c=1:numel(MeasuredData{1})
            if (strcmp(MeasuredData{1}(c).Name,'Root')==0) && (MeasuredData{1}(c).Total_Samples==0)
                t_file{d,1}=d;
                t_file{d,2}=MeasuredData{1}(c).Property(1).Value;
                d=d+1;
            end
        end
    else
        for c=1:numel(MeasuredData)
            if (strcmp(MeasuredData(c).Name,'Root')==0) && (MeasuredData(c).Total_Samples==0)
                t_file{d,1}=d;
                t_file{d,2}=MeasuredData(c).Property(1).Value;
                d=d+1;
            end
            
        end
    end
    
elseif numFiles>1
    d=1;
    for g=1:numFiles
        for h=1:numel(MeasuredData{g,1})
            if (strcmp(MeasuredData{g,1}(h).Name,'Root')==0) && (MeasuredData{g,1}(h).Total_Samples==0)
                t_file{d,1}=d;
                t_file{d,2}=MeasuredData{g,1}(h).Property(1).Value;
                d=d+1;
            end
        end
        
        
    end
end



end
