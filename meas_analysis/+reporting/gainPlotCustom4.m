function [p, rmse] = gainPlotCustom4(root, pattern, freq, location_str, loc_logic)
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

ft = fittype( 'A + B*(-sign(x-E)+1)*x*U^(x) + BB*(sign(x-E)+1)*x*UU^(x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf -Inf -Inf 0 0];
opts.StartPoint = [0.728661681678271 0.73784165379759 0.988417928784981 1.7 0 0.7069];

[p, gof] = fit(log10(R), G, ft, opts);
d_val= log10(logspace(log10(min(R)), log10(max(R)), 20));
g_val=feval(p,d_val);

% calculate the RMS error
% rmse = calcMse(R, G, p);
rmse = gof.rmse;

% plot the gains
h = figure();
semilogx(R,G, 'color', [0,0,0]+0.75, 'marker', 'o', 'linestyle' , 'none')

hold on
semilogx(10.^d_val,g_val,'b*-')
hold off       
legend('data','piecewise fit','Location','best')
xlabel('Distance, $d$ (m)','Interpreter','Latex');
ylabel('$10 \times log_{10}( \sum{|h(t)|^{2}} )$','Interpreter','Latex');
reporting.setCommonAxisProps();

end


function e = calcMse(R, G, p)
    G2=p(log10(R));
    e=sqrt(sum((G-G2).^2)/length(G));
end

