function stats2rfnestdp( obj, file_filter, stats_dir_path, emu_path )
%STATS2RFNESTDP Convert stats file reduced tap CIR to delay profile used by
% the RFNest emulator.
%
%   This function reads the MAT files in the STATS directory and produces
%   the a CIR for the RFNest emulator

if nargin < 3
    emu_path = '..\emu';
end
if nargin < 2
    stats_dir_path = '.';
end
if nargin < 1
    file_filter = '*.mat';
end

if ~exist(emu_path,'dir')
    mkdir(emu_path);
end

files = dir([stats_dir_path '\' file_filter]);
if isempty(files)
    warning 'no files found'
    return
end
NN = length(files);

for ii = 1:length(files)
    file_name = files(ii).name;
    stats_file_path = [stats_dir_path '\' file_name];
    stats = load(stats_file_path);
    stats = stats.stats;
    meta = stats.meta;
    
    % extract the sample rate
    Fs = meta.SampleRate_MHz_num*1e6;
    Ts = 1/Fs;
 
    % extract the average cir struct and write delay profile file
    fout_root = [emu_path '\' file_name(1:end-4)];
    fout_root = strrep(fout_root, '__channel_stats', '');
    for jj = 1:2
        acir_st = stats.avg_cir_st(jj);
        if ~isempty(acir_st.mag)
            ntaps = length(acir_st.mag);
            write_delay_profile(acir_st, ntaps, fout_root);
        end
    end
end

end

function write_delay_profile(acir_st, ntaps, fout_root)
    
    los = false;
    cir_class = cell2mat(acir_st.class);
    cir_t = acir_st.time;
    cir_mag = acir_st.mag;
    cir_angle = acir_st.angle;
    
    if strcmpi(cir_class,'los')
        los = true;
    end
    
    % open file for writing
    outfile_path = [fout_root '_' cir_class '.csv'];
    f = fopen(outfile_path,'w');
    if f == -1
        warning(['problem opending file for writing ' outfile_path])
        return
    end
    
    % line 1: channel type
    linev = ones(ntaps,1);
    if los
        linev(1) = 2;
    end
    line = sprintf('%d,', linev);
    line(end) = '';
    line = [line 13 10];
    fwrite(f, line);
    
    % line 2: Gain
    linev = floor(32767*cir_mag/sum(cir_mag));
    line = sprintf('%i,', linev);
    line(end) = '';
    line = [line 13 10];
    fwrite(f, line);
    
    % line 3: phase
    cir_angle = wrapTo2Pi(cir_angle)/(2*pi);
    linev = round(32767*cir_angle);
    line = sprintf('%i,', linev);
    line(end) = '';
    line = [line 13 10];
    fwrite(f, line);    
    
    % line 4: Doppler
    linev = zeros(ntaps,1);
    line = sprintf('%i,', linev);
    line(end) = '';
    line = [line 13 10];
    fwrite(f, line);    
    
    % line 5: Delay
    cir_t = 1e9*cir_t;
    linev = zeros(ntaps,1);
    linev(1) = cir_t(1); 
    for ii = 2:ntaps
        linev(ii) = cir_t(ii) - sum(linev(1:ii-1));
    end
    linev(1) = 16; %+16 is a work around for a HW bug in RFNest
    
    plot([cir_t(:) linev(:)],'+-'), legend('t','d')
    title(['delay ' cir_class ': ' fout_root],'interpreter','none')
    drawnow
    linev = round2(linev,4);
    line = sprintf('%.0f,', linev);
    line(end) = '';
    line = [line 13 10];
    fwrite(f, line);  
    
    fclose(f);
end

function z = round2(x,y)

if numel(y)>1
  error('Y must be scalar')
end
z = round(x/y)*y;

end





