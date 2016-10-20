function [p1, p2] = gainPlot(pattern,freq, location_str, loc_logic, infpt)
% Analyze complex impulse responses from measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

if nargin < 5
    infpt = 10;
end
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

p1 = [];
p2 = [];

%
% query the list of stats files
%
files = dir(pattern);
Nfiles = length(files);

% local variables
R = [];
G = [];

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
    
    % path gain
    R = [R; stats.path_gain_range_m(:)];
    G = [G; stats.path_gain_dB(:)];
    
end

if isempty(R)
    return
end

% compute the least squares fit of the gain data
G = G(~isnan(R));
R = R(~isnan(R));

% lower segment
Gfit1 = G(R<infpt);
Rfit1 = R(R<infpt);
p1 = polyfit(log10(Rfit1), Gfit1, 1);

% upper segment
Gfit2 = G(R>infpt);
Rfit2 = R(R>infpt);
p2 = polyfit(log10(Rfit2), Gfit2, 1);

% find draw points
x_intersect = fzero(@(x) polyval(p1-p2,x),[log10(min(Rfit1)), log10(max(Rfit2))]);
d_val1 = log10(logspace(log10(min(Rfit1)), x_intersect, 20));
g_val1 = polyval(p1, d_val1);
d_val2 = log10(logspace(x_intersect, log10(max(Rfit2)), 20));
g_val2 = polyval(p2, d_val2);

% plot the gains
semilogx(R,G, 'color', [0,0,0]+0.75, 'marker', '.', 'linestyle' , 'none')
hold on
semilogx(10.^d_val1,g_val1,'b+-')
semilogx(10.^d_val2,g_val2,'b+-')
hold off

end



