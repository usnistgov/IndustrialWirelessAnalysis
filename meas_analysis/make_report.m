% Produce Report on all measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

%% Summary Statistics
%summarize_stats();

%top_dirs = { 'test_report'};
%top_dirs = { 'AAplant', 'Boulder', 'GBurg'};
top_dirs = { 'Boulder'};

%% Plots of Measurements
for ii = 1:length(top_dirs)
    
    cd(top_dirs{ii})
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
        fprintf('Mat File: %s\n', meta.MatFile_str);
        fprintf('Frequency: %0.3f\n', meta.Frequency_GHz_num);
        fprintf('Location: %s\n', meta.Location_str);
        fprintf('RX Antenna: %s\n', meta.ReceiverAntenna_str);
        fprintf('RX Antenna Gain: %0.3f\n', meta.ReceiverAntennaGain_dBi_num);
        fprintf('TX Antenna: %s\n', meta.TransmitterAntenna_str);
        fprintf('TX Antenna Gain: %0.3f\n', meta.TransmitterAntennaGain_dBi_num);
        fprintf('TX Power, Watts: %0.3f\n', meta.TransmitterPower_watts_num);
        fprintf('PN Oversample Factor: %s\n', meta.PNOversample_num);
        fprintf('Sample rate, MHz: %0.3f\n', meta.SampleRate_MHz_num);
        fprintf('Frequency: %0.3f\n', meta.Frequency_GHz_num);

        % load plots
        fig_files = dir(['.\figs\' mat_root '*']);
        KK = length(fig_files);
        figs = [];
        for kk = 1:KK
            fig_fname = fig_files(kk).name;
            title_str = [];
            if ~isempty(strfind(fig_fname,'avgcir_los'))
                figs(3).title_str = 'Averge CIR for LOS Channels';
                figs(3).fig_fname = fig_fname;
            elseif ~isempty(strfind(fig_fname,'avgcir_nlos'))
                figs(4).title_str = 'Averge CIR for Non-LOS Channels';
                figs(4).fig_fname = fig_fname;
            elseif ~isempty(strfind(fig_fname,'__plvacq'))
                figs(1).title_str = 'Path Loss versus Acquisition Order';
                figs(1).fig_fname = fig_fname;
            elseif ~isempty(strfind(fig_fname,'__ds'))
                figs(6).title_str = 'Power Delay Spread';
                figs(6).fig_fname = fig_fname;
            elseif ~isempty(strfind(fig_fname,'__du'))
                figs(5).title_str = 'Power Delay';
                figs(5).fig_fname = fig_fname;
            elseif ~isempty(strfind(fig_fname,'__Kcdf'))
                figs(7).title_str = 'Rician K-factor';
                figs(7).fig_fname = fig_fname;
            elseif ~isempty(strfind(fig_fname,'__pl_dB'))
                figs(2).title_str = 'Path Loss versus Distance';
                figs(2).fig_fname = fig_fname;
            end
            
        end
        
        for kk = 1:length(figs)
            if ~isempty(figs(kk).fig_fname)
                fig_path = ['.\figs\' figs(kk).fig_fname];
                title_str = figs(kk).title_str;
                fig = openfig(fig_path);
                title(title_str);
                snapnow;
                close(fig)
            end
        end
    
    end
    
    % summary report for site
    
    cd('..\')

end





