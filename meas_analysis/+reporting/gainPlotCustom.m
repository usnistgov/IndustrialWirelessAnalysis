function [p1, p2] = gainPlotCustom(pattern,freq, location_str, loc_logic)
% Analyze complex impulse responses from measurements
% Author: Rick Candell, Mohamed Hany
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
        if strfind(stats_file_path,'Oats')
            continue;
        end
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

ft = fittype('A*exp(x)+B*x+C');
[p,S] = fit(log10(R), G, ft);
d_val= log10(logspace(log10(min(R)), log10(max(R)), 20));
g_val=feval(p,d_val);

% plot the gains
h = figure();
semilogx(R,G, 'color', [0,0,0]+0.75, 'marker', '.', 'linestyle' , 'none')

    hold on
    semilogx(10.^d_val,g_val,'b*-')
    hold off
    legend('Measured',sprintf('P: %0.3f log(x) + %0.3f x + %0.3f, Erms %.1f dB',p.A, p.B, p.C, 10*log10(S.rmse)),...
        'Location','southwest')
xlabel('Distance, d (m)','Interpreter','Latex')
ylabel('$10\ log_{10}(G)$','Interpreter','Latex')
ylim([-120 -40])
reporting.setCommonAxisProps()
reporting.setFigureForPrinting(h)

end


function e = calcMse(R, G, pcustom)
end

