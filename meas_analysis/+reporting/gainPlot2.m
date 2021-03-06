function [p1, p2, rmse] = gainPlot2(root, pattern, freq, location_str, loc_logic, infpt)
% Analyze complex impulse responses from measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

if nargin < 6
    infpt = 10;
end
if nargin < 5
    loc_logic = true;
end
if nargin < 4
    location_str = NaN;
end
if nargin < 3
    error 'frequency specification is required'
end
if nargin <2
    error 'pattern is required'
end
if nargin <1
    error 'root folder is required'
end

p1 = [];
p2 = [];
rmse = nan;

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
        
        stats_file_path = [root '\' files(fk).name];      
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
    R = [R; stats.path_gain_range_m(:)]; %#ok<AGROW>
    G = [G; stats.path_gain_dB(:)]; %#ok<AGROW>
    
end

if isempty(R)
    return
end

% compute the least squares fit of the gain data
G = G(~isnan(R));
R = R(~isnan(R));
% G=G(R>=10);
% R=R(R>=10);

% lower segment
Gfit1 = G(R<infpt);
Rfit1 = R(R<infpt);
[p1, S1] = polyfit(log10(Rfit1), Gfit1, 1);

% upper segment
Gfit2 = G(R>infpt);
Rfit2 = R(R>infpt);

% compute polynomial fit
[p2, S2] = polyfit(log10(Rfit2), Gfit2, 1);

% find draw points
legcnt = 3;
try
    x_intersect = fzero(@(x) polyval(p1-p2,x),[log10(min(Rfit1)), log10(max(Rfit2))]);
    d_val1 = log10(logspace(log10(min(Rfit1)), x_intersect, 20));
    [g_val1] = polyval(p1, d_val1, S1);
    d_val2 = log10(logspace(x_intersect, log10(max(Rfit2)), 20));
    [g_val2] = polyval(p2, d_val2, S2);
    rmse = calcMse(R, G, legcnt, x_intersect, p1, p2);
catch me
    legcnt = 2;
    d_val2 = log10(logspace(log10(infpt), log10(max(Rfit2)), 20));
    [g_val2] = polyval(p2, d_val2, S2);
    rmse = calcMse(R, G, legcnt, nan, p1, p2);
end

% plot the gains
h = figure();
semilogx(R,G, 'color', [0,0,0]+0.75, 'marker', '.', 'linestyle' , 'none')
if legcnt == 3
    hold on
    semilogx(10.^d_val1,g_val1,'b+-')
    semilogx(10.^d_val2,g_val2,'b*-')
    hold off    
    legend('Measured',sprintf('P1: %0.1fd + %0.1f, rmse %.1f dB',p1(2), p1(1),rmse), ...
        sprintf('P2: %0.1fd + %0.1f, rmse %.1f dB',p2(2), p2(1), rmse), ...
        'Location','southwest')
else
    hold on
    semilogx(10.^d_val2,g_val2,'b*-')
    hold off
    legend('Measured',sprintf('P2: %0.1f log10(d) + %0.1f, rmse %.1f dB',p2(2), p2(1), rmse),...
        'Location','southwest')
end
xlabel('Distance, d (m)','Interpreter','Latex')
ylabel('$10\ log_{10}(G)$','Interpreter','Latex')
ylim([-120 -40])
reporting.setCommonAxisProps()
reporting.setFigureForPrinting(h)

end


function e = calcMse(R, G, legcnt, x_intersect, p1, p2)
    if legcnt ==3
        G2_1=p1(1)*log10(R(R<x_intersect))+p1(2);
        G2_2=p2(1)*log10(R(R>=x_intersect))+p2(2);
        G2=[G2_1 G2_2];
        e=sqrt(sum((G-G2).^2)/length(G));
    else
        G2=p2(1)*log10(R)+p2(2);
        e=sqrt(sum((G-G2).^2)/length(G));
    end 
end

