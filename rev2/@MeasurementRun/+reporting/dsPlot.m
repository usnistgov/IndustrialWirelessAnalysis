function [Smean, Smax, Smin] = dsPlot(pattern,freq, location_str, loc_logic)
% Analyze complex impulse responses from measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

if nargin < 4
    loc_logic = true;
end
if nargin < 3
    location_str = NaN;
end
if nargin < 2
    error 'frequency specification is required'
end
if nargin <1
    error 'patter is required'
end

%
% query the list of stats files
%
files = dir(pattern);
Nfiles = length(files);

% local variables
S = [];

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
    
    disp(['Processing file: ' stats_file_path])
    
    % delay spread
    S = [S; stats.rms_delay_spread_sec(:)];
    
end

if isempty(S)
    return
end
S = S(~isnan(S));
S = S*1e9;

Smean = mean(S);
Smax = max(S);
Smin = min(S);
Smad = mad(S);
Sstd = std(S);

% plot the histogram of S   
X = S;
X(abs(X)>(Smean+100)) = [];
plotPDFCDF(X, 300, 'Delay Spread, S', 'ns');
setCommonAxisProps();

% 
% [ds_counts,ds_centers] = hist(S,100);
% ds_probs = cumsum(ds_counts/sum(ds_counts));
% ds_centers = ds_centers(ds_probs < 0.99);
% ds_probs = ds_probs(ds_probs < 0.99);
% yyaxis left, area(ds_centers, [0 diff(ds_probs)],'FaceAlpha',0.25)
% str = 'Pr. $$\hat{S}$$'; ylabel(str,'Interpreter','Latex');
% yyaxis right, plot(ds_centers,ds_probs)
% ylim([0 1])
% str = 'Pr. $$\hat{S} < S$$'; ylabel(str,'Interpreter','Latex');
% str = 'rms delay spread, $$S$$ (ns)';xlabel(str,'Interpreter','Latex')
% reporting.setCommonAxisProps();   
% reporting.setFigureForPrinting(gcf);

end



