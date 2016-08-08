function cmp_pl_poly( file_filter, stats_dir_path, fig_dir_path, png_dir_path )
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
    file_filter = '*.mat'
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

for ff = freqs(:)'
    run_names = {};
    d = linspace(1,300,100);
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
        if freq == ff
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
    run_names{end+1} = ['combined (' sprintf('p=%0.2fx + %0.2f',p_all) ')'];

    h = figure(gcf);
    semilogx(d, pv_arr, d, pv_all, 'k+-')
    text(10,max(pv_all)+10,run_names{end});
    xlabel('distance (m)')
    ylabel('Gain (dB)')

    savefig(h, [fig_dir_path '\zCombined_' num2str(ff) 'GHz__pathloss.fig']);
    print(h, [png_dir_path '\zCombined_' num2str(ff) 'GHz__pathloss.png'],'-dpng')
    close(h) 
    
end

end

