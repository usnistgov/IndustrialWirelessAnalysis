function Cstats = makeStatsSummary(out_path, pattern, freq, location_str, loc_logic)
% Analyze complex impulse responses from measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

if nargin < 5
    loc_logic = true;
end
if nargin < 4
    location_str = NaN;
end
if nargin < 3
    freq = NaN;
end
if nargin <2
    error 'pattern is required'
end
if nargin < 1
    error 'output file path is required'
end

%
% query the list of stats files
%
files = dir(pattern);
Nfiles = length(files);

% local variables
Cstats = {};

%
% Process each file in turn
%
for fk = 1:Nfiles
    
    % use test data or the real thing
    try 
        
        stats_file_path = files(fk).name;      
        stats = load(stats_file_path);
        stats = stats.stats;
        
    catch me
        warning('Problem reading mat file, trying again then skipping.');
        disp(me.message)         
        warning('Skipping file...');
        continue;
    end
    
    try
        meta = stats.meta;
    catch me
        warning('problem reading meta data')
        disp(me.message);
        continue;
    end
    
    % filter files
    if ~reporting.keepFile(meta, freq, location_str, loc_logic)
        continue;
    end
    
    disp(['Processing file: ' meta.MatFile_str])
    
    % save the stats for this measurement run
    RxPol = 'U';
    if strfind(meta.ReceiverAntenna_str,'V Pol')
        RxPol = 'V';
    elseif ~isempty(strfind(meta.ReceiverAntenna_str,'Cross Pol')) || ...
            ~isempty(strfind(meta.ReceiverAntenna_str,'X Pol')) || ...
            ~isempty(strfind(meta.ReceiverAntenna_str,'Long Pol'))
        RxPol = 'X';
    end
    TxPol = 'U';
    if ~isempty(strfind(meta.TransmitterAntenna_str,'V Pol'))
        TxPol = 'V';
    elseif ~isempty(strfind(meta.TransmitterAntenna_str,'Cross Pol')) || ...
            ~isempty(strfind(meta.TransmitterAntenna_str,'X Pol')) || ...
            ~isempty(strfind(meta.TransmitterAntenna_str,'Long Pol'))
        TxPol = 'X';
    end    
    Cstats(fk,:) = {   
        meta.MatFile_str, meta.Frequency_GHz_num, ...
        RxPol, meta.ReceiverAntennaGain_dBi_num, ...
        TxPol, meta.TransmitterAntennaGain_dBi_num, ...
        stats.path_gain_dB_poly(1), stats.path_gain_dB_poly(2), ...
        nanmean(stats.K), nanmin(stats.K), nanmax(stats.K) ...
        1e9*nanmean(stats.mean_delay_sec), 1e9*nanmin(stats.mean_delay_sec), 1e9*nanmax(stats.mean_delay_sec)...
        1e9*nanmean(stats.rms_delay_spread_sec), 1e9*nanmin(stats.rms_delay_spread_sec), 1e9*nanmax(stats.rms_delay_spread_sec) ...
    };
    
end

writeStatsToFile(Cstats, out_path);

end


function writeStatsToFile(stats_arr, file_path)
    if isempty(stats_arr)
        return
    end
    M = cell2table(stats_arr);
    VariableNames = ...
        {'Run', 'Freq', 'RX_Ant_Type', 'RX_Ant_Gain', 'TX_Ant_Type', 'TX_Ant_Gain', ...
        'Path_Gain_Poly_Slope', 'Path_Gain_Poly_YInt', 'Mean_K', 'Min_K', 'Max_K', ...
        'Mean_Delay_ns', 'Min_Delay_ns', 'Max_Delay_ns', ...
        'Mean_Delay_Spread_ns', 'Min_Delay_Spread_ns', 'Max_Delay_Spread_ns'};
    M.Properties.VariableNames = VariableNames;
    
    M0 = [];
    if exist(file_path,'file')
        M0 = readtable(file_path);
    end
    if ~isempty(M0)
        M = [M0;M];
    end
    M.Properties.VariableNames = VariableNames;
    writetable(M, file_path, 'Delimiter', '\t');
    
end


