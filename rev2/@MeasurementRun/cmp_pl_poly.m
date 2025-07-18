function cmp_pl_poly(obj, file_filter, stats_dir_path, fig_dir_path, png_dir_path )
%CMP_PL_POLY Compare pathloss across measurement runs
%   This function reads the MAT files in the STATS directory and compares
%   the path loss polynomial fits

if nargin < 4
    png_dir_path = '..\png';
end
if nargin < 3
    fig_dir_path = '..\figs';
end
if nargin < 2
    stats_dir_path = '.';
end
if nargin < 1
    file_filter = '*_stats.mat';
end

files = dir([stats_dir_path '\' file_filter]);
if isempty(files)
    warning 'no files found'
    return
end
NN = length(files);

freqs = nan(NN,1);
for jj = 1:length(files)
    file_name = files(jj).name;
    stats = load([stats_dir_path '\' file_name]);
    meta = stats.stats.meta;
    freqs(jj) = meta.Frequency_GHz_num;
end
freqs = unique(freqs);

% frequency by frequency
for ff = freqs(:)'
    plotByFreq(ff, files, stats_dir_path, fig_dir_path, png_dir_path);
end

end

function plotByFreq(ff, files, stats_dir_path, fig_dir_path, png_dir_path, title_str, show_traces)

    if nargin < 7
        show_traces = true;
    end

    if nargin < 6
        title_str = [];
    end

    run_names = {};
    %d = linspace(1,150,100);
    d = log10(logspace(log10(30),log10(200),10));
    polys_str = {};
    pv_arr = [];
    d_all = [];
    pv_all = [];
    jj = 0;
    for ii = 1:length(files)
        file_name = files(ii).name;
        stats = load([stats_dir_path '\' file_name]);
        stats = stats.stats;
        meta = stats.meta;
        freq = meta.Frequency_GHz_num;
        if any(freq == ff)
            jj = jj + 1;
            run_names{jj} = meta.MatFile_str; %#ok<*AGROW>
            disp(run_names{jj})
            p = stats.path_gain_dB_poly;
            pv = polyval(p,d);
            pv_arr(:,jj) = pv(:);
            polys_str{jj} = sprintf('p=%0.2fx + %0.2f',p);
            d_all = [d_all(:); d(:)];
            pv_all = [pv_all(:); pv(:)];
        end
    end

    p_all = polyfit(d_all, pv_all, 1);
    pv_all = polyval(p_all, d);
    run_names{end+1} = sprintf('p=%0.2fx + %0.2f',p_all);

    h = figure(gcf);
    if show_traces
        semilogx(10.^d, pv_arr, 'color', [0,0,0]+0.7);
    end
    hold on
    % plot the traces
    semilogx(10.^d, pv_all, 'k+-')
    text(2,max(pv_all)-20,run_names{end});
    
    % frii as reference
    if length(ff) == 1
        d_frii = logspace(min(d), max(d), 20);
        c = physconst('LightSpeed');
        frii_fspl_dB = 10*log10(d_frii.^2) + 20*log10(ff) + 20*log10(1e9) + 20*log10(4*pi/c);
        semilogx(d_frii, -frii_fspl_dB, 'b-')
        text(d_frii(floor(length(d)/4)),-frii_fspl_dB(floor(length(d_frii)/4))+5,'FSPL','Color','blue');
    end
    
    hold off
%     if isempty(title_str)
%         title(['Channel Gain for ' sprintf('%0.3f GHz ', ff)])
%     else
%         title(title_str)
%     end
    xlabel('distance (m)')
    ylabel('Gain (dB)')

    setCommonAxisProps();
    if isempty(title_str)
        savefig(h, [fig_dir_path '\zCombined_' num2str(ff) 'GHz__pathloss.fig']);
        setFigureForPrinting();
        print(h, [png_dir_path '\zCombined_' num2str(ff) 'GHz__pathloss.png'],'-dpng')
    else
        savefig(h, [fig_dir_path '\zCombined_allfreqs_pathloss.fig']);
        setFigureForPrinting();
        print(h, [png_dir_path '\zCombined_allfreqs_pathloss.png'],'-dpng')        
    end
%     close(h) 

end

function setFigureForPrinting()
    width=3; height=3;
    set(gcf,'InvertHardcopy','on');
    set(gcf,'PaperUnits', 'inches');
    papersize = get(gcf, 'PaperSize');
    left = (papersize(1)- width)/2;
    bottom = (papersize(2)- height)/2;
    myfiguresize = [left, bottom, width, height];
    set(gcf,'PaperPosition', myfiguresize);
end

function setCommonAxisProps()

    alw = 0.75;    % AxesLineWidth
    fsz = 10;      % Fontsize
    lw = 1.5;      % LineWidth
    msz = 3.5;     % MarkerSize
    
%    grid on
    set(gca,'XGrid','on')
    set(gca,'XMinorGrid','off')
    set(gca,'YGrid','on')
    set(gca,'YMinorGrid','off')
    set(gca,'GridAlpha',0.5)
    set(gca,'MinorGridAlpha',0.4)
    set(gca,'Fontsize',fsz)
    set(gca,'LineWidth',alw);
    set(gca,'FontName','TimesRoman')
    
    % set the line properties
    hline = get(gca,'Children');
    for h = hline(:)'
        if isa(h,'matlab.graphics.chart.primitive.Line')
            h.LineWidth = lw;
            h.MarkerSize = msz;
        elseif isa(h,'matlab.graphics.primitive.Text')
            set(h,'FontName','TimesRoman')
        end
    end
end
