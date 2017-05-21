function makeLatexFigReport(fname)
% Produce Report on all measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

if nargin < 1
    fname = '../pub/deltailedresults.tex';
end

disp(['output going to: ' fname])

%% Summary Statistics
% summarize_stats();

% output stream
f = fopen(fname,'w');  % standard output
fprintf('%automatically generated from amtlab on %s', datestr(datetime('now')));
stats_files = dir('.\*_stats.mat');
NN = length(stats_files);
    
for jj = 1:NN

    mat_fname = stats_files(jj).name;
    mat_root = strrep(mat_fname,'_pp__channel_stats.mat','');

    % load the stats file
    stats_fname = stats_files(jj).name;
    stats = load(stats_fname);
    stats = stats.stats;
    meta = stats.meta;
    
    % add section title
    addSubSubSection(f, meta.MatFile_str);

    % summary tables
    addParametersTable(f, meta);
    addMetricsTable(f, stats, true);
    addMetricsTable(f, stats, false);

    % load PNG plots
    png_files = dir(['..\png\' mat_root '*']);
    KK = length(png_files);
    pngs = [];
    for kk = 1:KK
        png_fname = png_files(kk).name;
        title_str = [];
        mat_root_str = strrep(mat_root, '_', ' ');
        if ~isempty(strfind(png_fname,'avgcir_los'))
            %pngs(3).title_str = sprintf('Average CIR for LOS Channels (%s)', mat_root_str);
            %pngs(3).png_fname = png_fname;
        elseif ~isempty(strfind(png_fname,'avgcir_nlos'))
            %pngs(4).title_str = sprintf('Average CIR for Non-LOS Channels (%s)', mat_root_str);
            %pngs(4).png_fname = png_fname;
        elseif ~isempty(strfind(png_fname,'__plvacq'))
            pngs(1).title_str = sprintf('Path Gain versus Acquisition Order');
            pngs(1).png_fname = png_fname;
        elseif ~isempty(strfind(png_fname,'__ds'))
            pngs(6).title_str = sprintf('Power Delay Spread');
            pngs(6).png_fname = png_fname;
        elseif ~isempty(strfind(png_fname,'__du'))
            pngs(5).title_str = sprintf('Power Delay');
            pngs(5).png_fname = png_fname;
        elseif ~isempty(strfind(png_fname,'__Kcdf'))
            pngs(7).title_str = sprintf('Rician K-factor');
            pngs(7).png_fname = png_fname;
        elseif ~isempty(strfind(png_fname,'__pl'))
            extra = sprintf('The effective path loss exponent of this measurement dataset is %0.2f.', abs(stats.path_gain_dB_poly(1))/10);
            pngs(2).title_str = sprintf('Path Loss versus Distance');
            pngs(2).png_fname = png_fname;
        end

    end
     
    % Figure setion header
    %addParagraphSection(f, 'Channel Estimation Analysis');
    
    % add each figure separately in the document
    counter = 0;
    fprintf(f, '\\begin{figure}[H]\n');
    fprintf(f, '\\centering\n');
    for kk = 1:length(pngs)
        if ~isempty(pngs(kk).png_fname)
            png_path = ['../png/' pngs(kk).png_fname];
            pub_png_path = ['../pub/figs/' pngs(kk).png_fname];
            pub_png_path = strrep(pub_png_path, '_', '');
            pub_png_path = strrep(pub_png_path, ' ', '');
            tex_png_path = ['figs/' pngs(kk).png_fname(1:end-4)];
            tex_png_path = strrep(tex_png_path, '_', '');
            tex_png_path = strrep(tex_png_path, ' ', '');
            copyfile(png_path, pub_png_path);
            s_cap = strrep(pngs(kk).title_str, '_', ' ');
            addPngFigure(f, tex_png_path, s_cap);    
            if counter
                fprintf(f, '\\ \\ \\ \\ \n');
                counter = 0;
            else
                counter = counter + 1;
            end
        end
    end
    capstr = sprintf('Channel characterizations for %s', strrep(mat_fname(1:end-4),'_',' '));
    fprintf(f, '\\caption{%s}\n', capstr);
    fprintf(f, '\\end{figure}\n\n');
    
    addNewPage(f);

end

fclose(f);

end % function

function addSubSection(f, name)
    fprintf(f, '\\subsection{%s}\n', strrep(name,'_',' ') );
end

function addSubSubSection(f, name)
    fprintf(f, '\\subsubsection{%s}\n', strrep(name,'_',' ') );
end

function addParagraphSection(f, name)
    fprintf(f, '\\paragraph{%s}\n', strrep(name,'_',' ') );
end

function addNewPage(f)
    fprintf(f, '\\newpage\n');
end

function addPngFigure(f, png_path, s_cap)
%     fprintf(f, '\\begin{figure}[H]\n');
    fprintf(f, '\n\\subfloat[%s]\n{\n', s_cap);
%     fprintf(f, '\\centering\n');
%     fprintf(f, '\\captionsetup{width=3in}\n');
    png_path = strrep(png_path,'\','/');
    png_path = strrep(png_path,'_','_');
    fprintf(f, '   \\includegraphics[width=2.5in]{%s} ', png_path);
    fprintf(f,'    \\label{fig:detailed:%s}\n',strrep(png_path,'/',':') );    
%     fprintf(f, '\\caption{%s}\n', s_cap);
    fprintf(f, '}');
%     fprintf(f, '\\end{figure}\n');
end

function addParametersTable(f, meta)
    
        fprintf(f, '\n');
        fprintf(f, '\\begin{table}[H]\n');
        fprintf(f, '\\renewcommand{\\arraystretch}{1.3}\n');
        fprintf(f, '\\centering\n');
        fprintf(f, '\\caption{Measurement Parameters}\n');
        %fprintf(f, '\\label{results:params_%s}\n', meta.MatFile_str);
        
        fprintf(f, '\\begin{tabular}{l|c}\n');
        fprintf(f, '\\hline\n');
        fprintf(f, '\\textbf{Parameter}&\\textbf{Value}\\\\\\hline\\hline\n');
        fprintf(f, 'Mat File & %s\\\\\n', strrep(meta.MatFile_str,'_','\_'));
        fprintf(f, 'Frequency (GHz) & %0.3f\\\\\n', meta.Frequency_GHz_num);
        fprintf(f, 'Location & %s\\\\\n', meta.Location_str);
        fprintf(f, 'RX Antenna & %s\\\\\n', meta.ReceiverAntenna_str);
        fprintf(f, 'RX Antenna Gain & %0.3f\\\\\n', meta.ReceiverAntennaGain_dBi_num);
        fprintf(f, 'TX Antenna & %s\\\\\n', meta.TransmitterAntenna_str);
        fprintf(f, 'TX Antenna Gain & %0.3f\\\\\n', meta.TransmitterAntennaGain_dBi_num);
        fprintf(f, 'TX Power, Watts & %0.3f\\\\\n', meta.TransmitterPower_watts_num);
        fprintf(f, 'PN Oversample Factor & %0.3f\\\\\n', meta.PNOversample_num);
        fprintf(f, 'Sample rate, MHz & %0.3f\\\\\n', meta.SampleRate_MHz_num);
        fprintf(f, '\\end{tabular}\n');
        fprintf(f, '\\end{table}\n\n');
end

function addMetricsTable(f, stats, delol)

    if nargin < 3
        delol = false;
    end

    K = stats.K;
    K       = K(~isnan(K));
    if delol
       K = deleteoutliers(K);
    end    
    meanK   = mean(K);
    medianK = median(K);
    minK    = min(K);
    maxK    = max(K);
    madK    = mad(K);    
    stdK    = std(K);
    
    Klos = stats.K(stats.los==1);
    Klos       = Klos(~isnan(Klos));
    if delol
       Klos = deleteoutliers(Klos);
    end        
    meanKlos   = mean(Klos);
    medianKlos = median(Klos);
    minKlos    = min(Klos);
    maxKlos    = max(Klos);
    madKlos    = mad(Klos);    
    stdKlos    = std(Klos);    
    
    Tau = 1e9*stats.mean_delay_sec;
    Tau     = Tau(~isnan(Tau));
    if delol
       Tau = deleteoutliers(Tau);
    end      
    meanTau = mean(Tau);
    medianTau = median(Tau);
    minTau  = min(Tau);
    maxTau  = max(Tau);
    madTau  = mad(Tau);
    stdTau  = std(Tau);
    
    S = 1e9*stats.rms_delay_spread_sec;
    S       = S(~isnan(S));
    if delol
       S = deleteoutliers(S);
    end      
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
    if delol
        fprintf(f, '\\caption{Statistics of Channel Estimates.  Outliers are removed using a significance test of 0.95.}\n');
    else
        fprintf(f, '\\caption{Statistics of Channel Estimates.  Outliers are not removed.  N indicates the total number of samples in the population.  Outliers are provided at 10 and 20 times $\\sigma$. }\n');
    end
    
    fprintf(f, '\\begin{tabular}{l||c|c|c|c||c|c||l}\n');
    fprintf(f, '\\hline\n');
    fprintf(f, 'Estimate        & min$(\\cdot)$   & median$(\\cdot)$ & mean$(\\cdot)$     & max$(\\cdot)$   & std$(\\cdot)$   & mad$(\\cdot)$   & outlier info \\\\\\hline\\hline\n');
    
    % K
    fprintf(f, '$K$ (dB),~NLOS  & %0.1f & %0.1f  & %0.1f    & %0.1f & %0.1f & %0.1f & N: %i   \\\\\n', ...
                                minK,   medianK,    meanK,  maxK,   stdK,   madK,   length(K));
	addOutlierEstimates(f,K,meanK,stdK);
    fprintf(f, '\\hline\n');
    
    % K los
    fprintf(f, '$K$ (dB),~LOS  & %0.1f & %0.1f  & %0.1f    & %0.1f & %0.1f & %0.1f & N: %i   \\\\\n', ...
                                minKlos,   medianKlos,    meanKlos,  maxKlos,   stdKlos,   madKlos,   length(Klos));
	addOutlierEstimates(f,Klos,meanKlos,stdKlos);
    fprintf(f, '\\hline\n');   
    
    %Tau
    fprintf(f, '$\\tau$ (ns)    & %0.1f & %0.1f  & %0.1f    & %0.1f & %0.1f & %0.1f & N: %i   \\\\\n', ...
                                minTau,   medianTau,    meanTau,  maxTau,   stdTau,   madTau,   length(Tau));
	addOutlierEstimates(f,Tau,meanTau,stdTau);
    fprintf(f, '\\hline\n');
    
    % S
    fprintf(f, '$S$ (ns)    & %0.1f & %0.1f  & %0.1f    & %0.1f & %0.1f & %0.1f & N: %i   \\\\\n', ...
                                minS,   medianS,    meanS,  maxS,   stdS,   madS,   length(S));
	addOutlierEstimates(f,S,meanS,stdS);
    fprintf(f, '\\hline\n');
   
    fprintf(f, '\\end{tabular}\n');
    fprintf(f, '\\end{table}\n\n');    

end

function addOutlierEstimates(f,x,u,s)
	for nsig = 10:10:20
        nout = length(x(x>(u+nsig*s)));
        if nout
            fprintf(f, ' & & & & & & & $%i\\sigma (%0.1f): %i$ \\\\\n', ...
            nsig, u+nsig*s,  nout);
        else
            break;
        end
    end
end
