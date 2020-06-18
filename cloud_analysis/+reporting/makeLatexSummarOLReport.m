function [Kmean, Kmax, Kmin] = makeLatexSummarOLReport(ofname, pattern,freq, location_str, loc_logic)
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
    error 'frequency specification is required'
end
if nargin <2
    error 'pattern is required'
end
if nargin <1
    error 'output file is required'
end

%
% query the list of stats files
%
files = dir(pattern);
Nfiles = length(files);

% local variables
K = [];

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
    
    % delay spread
    K = [K; stats.K(:)];
    
end

if isempty(K)
    return
end
K = K(~isnan(K));

Kmean = mean(K);
Kmax = max(K);
Kmin = min(K);

tbl = makeSummaryOutlierTable('summwithol.tex', stats);
matrix2latex(tbl(1:N2,:), ofname, 'columnLabels', ...
    {'Run','$f$ (GHz)','$\gamma_2$', ...
        '$\mathbf{E}(\tau)$ (ns)','$\mathbf{E}(S)$ (ns)','$\mathbf{E}(K)$ (dB)',...
        'Max K (dB)'}, 'format', '%-.1f', 'alignment', 'c');

end


function makeSummaryOutlierTable(stats)

    %K = deleteoutliers(stats.K);
    K = stats.K;
    K       = K(~isnan(K));
    meanK   = mean(K);
    medianK = median(K);
    minK    = min(K);
    maxK    = max(K);
    madK    = mad(K);    
    stdK    = std(K);
    
    Klos = stats.K(stats.los==1);
    Klos       = Klos(~isnan(Klos));
    meanKlos   = mean(Klos);
    medianKlos = median(Klos);
    minKlos    = min(Klos);
    maxKlos    = max(Klos);
    madKlos    = mad(Klos);    
    stdKlos    = std(Klos);    
    
    %Tau = deleteoutliers(1e9*stats.mean_delay_sec);
    Tau = 1e9*stats.mean_delay_sec;
    Tau     = Tau(~isnan(Tau));
    meanTau = mean(Tau);
    medianTau = median(Tau);
    minTau  = min(Tau);
    maxTau  = max(Tau);
    madTau  = mad(Tau);
    stdTau  = std(Tau);
    
    %S = deleteoutliers(1e9*stats.rms_delay_spread_sec);
    S = 1e9*stats.rms_delay_spread_sec;
    S       = S(~isnan(S));
    meanS   = mean(S);
    medianS = median(S);
    minS    = min(S);
    maxS    = max(S);
    madS    = mad(S);
    stdS    = std(S);
    
    fprintf(f, '\n');
    fprintf(f, '\\begin{table}[H]\n');
    fprintf(f, '\\renewcommand{\\arraystretch}{1.3}\n');
    fprintf(f, '\\centering\n');
    fprintf(f, '\\caption{Statistics of Channel Estimates.}\n');
    
    fprintf(f, '\\begin{tabular}{l||c|c|c|c||c|c||l}\n');
    fprintf(f, '\\hline\n');
    fprintf(f, 'Estimate        & min$(\\cdot)$   & median$(\\cdot)$ & mean$(\\cdot)$     & max$(\\cdot)$   & std$(\\cdot)$   & mad$(\\cdot)$   & outliers \\\\\\hline\\hline\n');
    
    % K
    fprintf(f, '$K$ (dB),~NLOS  & %0.1f & %0.1f  & %0.1f    & %0.1f & %0.1f & %0.1f & numel: %i   \\\\\n', ...
                                minK,   medianK,    meanK,  maxK,   stdK,   madK,   length(K));
	addOutlierEstimates(f,K,meanK,stdK);
    fprintf(f, '\\hline\\hline\n');
    
    % K los
    fprintf(f, '$K$ (dB),~LOS  & %0.1f & %0.1f  & %0.1f    & %0.1f & %0.1f & %0.1f & numel: %i   \\\\\n', ...
                                minKlos,   medianKlos,    meanKlos,  maxKlos,   stdKlos,   madKlos,   length(Klos));
	addOutlierEstimates(f,Klos,meanKlos,stdKlos);
    fprintf(f, '\\hline\\hline\n');   
    
    %Tau
    fprintf(f, '$\\tau$ (ns)    & %0.1f & %0.1f  & %0.1f    & %0.1f & %0.1f & %0.1f & numel: %i   \\\\\n', ...
                                minTau,   medianTau,    meanTau,  maxTau,   stdTau,   madTau,   length(Tau));
	addOutlierEstimates(f,Tau,meanTau,stdTau);
    fprintf(f, '\\hline\\hline\n');
    
    % S
    fprintf(f, '$S$ (ns)    & %0.1f & %0.1f  & %0.1f    & %0.1f & %0.1f & %0.1f & numel: %i   \\\\\n', ...
                                minS,   medianS,    meanS,  maxS,   stdS,   madS,   length(S));
	addOutlierEstimates(f,S,meanS,stdS);
    fprintf(f, '\\hline\\hline\n');
   
    fprintf(f, '\\end{tabular}\n');
    fprintf(f, '\\end{table}\n\n');    

end

function addOutlierEstimates(f,x,u,s)
	for nsig = 3:3:12
        nout = length(x(x>(u+nsig*s)));
        if nout
            fprintf(f, ' & & & & & & & $%i\\sigma (%0.1f): %i$ \\\\\n', ...
            nsig, u+nsig*s,  nout);
        else
            break;
        end
    end
end