function make_latex_report(fname, doc_root)
% Produce Report on all measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

if nargin < 1
    fname = 'E:\FACTORY_MEASUREMENTS2\active_meas\pub\test_report.tex';
end

%% Summary Statistics
% summarize_stats();

top_dirs = { 'test_report'};
%top_dirs = { 'AAplant', 'Boulder', 'GBurg'};
%top_dirs = { 'Boulder'};

% output stream
f = fopen(fname,'w');  % standard output

%% Plots of Measurements
for ii = 1:length(top_dirs)
    
    this_dir = top_dirs{ii};
    cd(this_dir);
    mat_files = dir('.\*_pp.mat');
    NN = length(mat_files);
    
    for jj = 1:NN
        
        mat_fname = mat_files(jj).name;
        mat_root = strrep(mat_fname,'_pp.mat','');
        
        % load the stats file
        stats_fname = ['.\stats\' mat_root '_pp__channel_stats.mat'];
        if ~exist(stats_fname,'file')
            continue;
        end
        stats = load(stats_fname);
        meta = stats.stats.meta;
        
        % summary table
        makeParametersTable(f, meta);

        % load PNG plots
        png_files = dir(['.\png\' mat_root '*']);
        KK = length(png_files);
        pngs = [];
        for kk = 1:KK
            png_fname = png_files(kk).name;
            title_str = [];
            if ~isempty(strfind(png_fname,'avgcir_los'))
                pngs(3).title_str = 'Averge CIR for LOS Channels';
                pngs(3).png_fname = png_fname;
            elseif ~isempty(strfind(png_fname,'avgcir_nlos'))
                pngs(4).title_str = 'Averge CIR for Non-LOS Channels';
                pngs(4).png_fname = png_fname;
            elseif ~isempty(strfind(png_fname,'__plvacq'))
                pngs(1).title_str = 'Path Loss versus Acquisition Order';
                pngs(1).png_fname = png_fname;
            elseif ~isempty(strfind(png_fname,'__ds'))
                pngs(6).title_str = 'Power Delay Spread';
                pngs(6).png_fname = png_fname;
            elseif ~isempty(strfind(png_fname,'__du'))
                pngs(5).title_str = 'Power Delay';
                pngs(5).png_fname = png_fname;
            elseif ~isempty(strfind(png_fname,'__Kcdf'))
                pngs(7).title_str = 'Rician K-factor';
                pngs(7).png_fname = png_fname;
            elseif ~isempty(strfind(png_fname,'__pl_dB'))
                pngs(2).title_str = 'Path Loss versus Distance';
                pngs(2).png_fname = png_fname;
            end
            
        end
        
        for kk = 1:length(pngs)
            if ~isempty(pngs(kk).png_fname)
                png_path = ['png\' pngs(kk).png_fname];
                s_cap = pngs(kk).title_str;
                png_file2 = strrep(pngs(kk).png_fname,'_','');
                png_path2 = ['png\' png_file2];
                copyfile(png_path, ['../pub/png/' png_file2]);
                addPngFigure(f, png_path2, s_cap);            
%                 png = openpng(png_path);
%                 title(title_str);
%                 snapnow;
%                 close(png)
            end
        end
    
    end
    
    % summary report for site
    
    cd('..\')
    
end

fclose(f);

end % function

function addPngFigure(f, png_path, s_cap)
    fprintf(f, '\\begin{figure*}[!t]\n');
    fprintf(f, '\\centering\n');
    png_path = strrep(png_path,'_','_');
    png_path = strrep(png_path,'\','/');
    fprintf(f, '\\includegraphics[width=0.75\\linewidth]{%s}\n', png_path);
    fprintf(f, '\\caption{%s}\n', s_cap);
    fprintf(f, '\\end{figure*}\n\n');
end

function makeParametersTable(f, meta)
    
        fprintf(f, '\n');
        fprintf(f, '\\begin{table}[!t]\n');
        fprintf(f, '\\renewcommand{\\arraystretch}{1.3}\n');
        fprintf(f, '\\centering\n');
        fprintf(f, '\\caption{Measurement Parameters}\n');
        fprintf(f, '\\label{results:params_%s}\n', meta.MatFile_str);
        
        fprintf(f, '\\begin{tabular}{l|c}\n');
        fprintf(f, '\\hline\n');
        fprintf(f, '\\textbf{Parameter}&\\textbf{Value}\\\\\\hline\\hline\n');
        fprintf(f, 'Mat File & %s\\\\\n', strrep(meta.MatFile_str,'_','\_'));
        fprintf(f, 'Frequency & %0.3f\\\\\n', meta.Frequency_GHz_num);
        fprintf(f, 'Location & %s\\\\\n', meta.Location_str);
        fprintf(f, 'RX Antenna & %s\\\\\n', meta.ReceiverAntenna_str);
        fprintf(f, 'RX Antenna Gain & %0.3f\\\\\n', meta.ReceiverAntennaGain_dBi_num);
        fprintf(f, 'TX Antenna & %s\\\\\n', meta.TransmitterAntenna_str);
        fprintf(f, 'TX Antenna Gain & %0.3f\\\\\n', meta.TransmitterAntennaGain_dBi_num);
        fprintf(f, 'TX Power, Watts & %0.3f\\\\\n', meta.TransmitterPower_watts_num);
        fprintf(f, 'PN Oversample Factor & %0.3f\\\\\n', meta.PNOversample_num);
        fprintf(f, 'Sample rate, MHz & %0.3f\\\\\n', meta.SampleRate_MHz_num);
        fprintf(f, 'Frequency & %0.3f\\\\\n', meta.Frequency_GHz_num);
        fprintf(f, '\\end{tabular}\n');
        fprintf(f, '\\end{table}\n\n');
end



