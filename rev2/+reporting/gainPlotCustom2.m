function [p, rmse] = gainPlotCustom2(root, pattern, freq, location_str, loc_logic)
% Analyze complex impulse responses from measurements
% Author: Rick Candell, Mohamed Hany
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

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
    error 'patter is required'
end
if nargin <1
    error 'root folder is required'
end

p = [];
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
    R = [R; stats.path_gain_range_m(:)];
    G = [G; stats.path_gain_dB(:)];
    
end

if isempty(R)
    return
end

% compute the least squares fit of the gain data
G = G(~isnan(R));
R = R(~isnan(R));
G=G(R>=10);
R=R(R>=10);

ft = fittype('A*F^(x)+B*x+C');
ftvars = coeffnames(ft);
fo = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-50,-100, -100, 0],...
        'Upper',[  0, 20, 100,   20], ...
        'StartPoint',[-4 10 -50 exp(1)]);
p = fit(log10(R), G, ft, fo);
d_val= log10(logspace(log10(min(R)), log10(max(R)), 20));
g_val=feval(p,d_val);

% calculate the RMS error
rmse = calcMse(R, G, p);

% plot the gains
h = figure();
semilogx(R,G, 'color', [0,0,0]+0.75, 'marker', '.', 'linestyle' , 'none')

    hold on
    semilogx(10.^d_val,g_val,'b*-')
    hold off
    legend('Measured',...
        sprintf('P: %0.1f*%0.1f^{log10(d)} + %0.1f*log10(d) %+0.1f, rmse %.1f dB',p.A, p.F, p.B, p.C, rmse), ...
        'Location','southwest', 'Interpreter','Latex');
xlabel('Distance, d (m)','Interpreter','Latex');
ylabel('$10\ log_{10}(G)$','Interpreter','Latex');
ylim([-120 -40]);
reporting.setCommonAxisProps();
reporting.setFigureForPrinting(h);

end


function e = calcMse(R, G, pcustom)
    G2=pcustom.A*pcustom.F.^(log10(R))+pcustom.B*log10(R)+pcustom.C;
    e=sqrt(sum((G-G2).^2)/length(G));
end

