function [ h ] = plot_hilbert_analysis_eeg(params, summary, save_to_file)
%PLOT_HILBERT_SUMMARY_EEG Plots Hilbert summary statistics. Requires that
% analyze_hilbert_eeg.m has already been run.

if nargin < 3
    save2file = true;
end

N_subj = length(summary.subjects);
N_channels = length(summary.channels);
N_bands = height(summary.bands);

% 1. Epochs
%    Does Hilbert envelope magnitude differ between epochs? 

show_plots = params.eeg.hilbert.plot.show_plots || ~save_to_file;
plotly_layouts = load(params.general.plots.plotly_layouts_file);

outdir = sprintf('%s/figures', params.io.results_dir);
if ~exist(outdir, 'dir')
   mkdir(outdir); 
end

%    1.1. Passing vs. baseline

%    1.1.1 Boxplots - average over channels
h = plotlyfig('Visible','off');
h.data = cell(3,1);
Y = zeros(N_bands*N_subj,2);
X = cell(N_bands*N_subj,1);
G = [];
idx = 0;

% Boxplots groups by frequency band
for bb = 1 : N_bands
    band = summary.bands.Band{bb};
    result = summary.stats.ttest.grand_mean.(band).baseline_overtakes;
    
    Y(idx+1:idx+N_subj,:) = result.data;
    X(idx+1:idx+N_subj) = {band};

    G = result.labels;
    idx = idx + N_subj;
end

for i = 1 : 2
    h.data(i) = {struct( ...
                    'y', Y(:,i), ...
                    'x', {X}, ...
                    'name', G(i), ...
                    'type', 'box' ...
                    )};
end

% Significance asterisks
Y = ones(N_bands,1)*0.15;
X = {};
T = {};

for bb = 1 : N_bands
    band = summary.bands.Band{bb};
    result = summary.stats.ttest.grand_mean.(band).baseline_overtakes;
    if result.sig_fdr
        X = [X {band}];
        T = [T {'-*-'}];
    end
end

h.data(3) = {struct( ...
                'y', Y, ...
                'x', {X}, ...
                'name', 'p<0.05', ...
                'type', 'scatter', ...
                'mode', 'text', ...
                'textposition', 'top', ...
                'text', {T}, ...
                'showlegend', false, ...
                'textfont', struct(...
                      'family', 'Courier New, monospace', ...
                      'size', 40, ...
                      'color', 'black')...
                )};


h.layout = plotly_layouts.layouts.hilbert.boxplot;
h.layout.title = 'Hilbert envelopes: Baseline vs. Passing';

h.PlotOptions.SaveFolder = outdir;
h.PlotOptions.FileName = 'eeg_hilbert_bl_overtake_boxplots';

% try
    plotlyoffline(h);
% catch ex
%     fprintf('\nAttempting to fetch plotlyoffline dependencies...');
%     getplotlyoffline(params.general.plotly_url);
%     plotlyoffline(h);
% end

if show_plots
    web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
end


%    1.1.3 Topoplots - Baseline vs. Overttake - T-stat by channel

% Load layout
cfg = [];
cfg.layout = params.eeg.layout;
evalc('layout=ft_prepare_layout(cfg);');
cfg.layout = layout;
cfg.colorbar = 'yes';
cfg.marker = 'no';
cfg.zlim = [-5 5]; % params.eeg.hilbert.plot.topo.baseline_overtake.clim;
cfg.colorbar = 'yes';
cfg.colorbartext = 'T-stat';
cfg.colormap = '-RdBu'; % params.eeg.hilbert.plot.topo.baseline_overtake.colormap;

h = figure('Visible','off');
h.Color = 'w';

h.Position([3 4]) = [320*N_bands 250];
t = tiledlayout("horizontal");
title(t, 'Hilbert transform (Passing-Baseline)');

for bb = 1 : N_bands
    nexttile;

    % Prepare stats
    band = summary.bands.Band{bb};
    result = summary.stats.ttest.channels.(band).baseline_overtakes;
    
    % Plot t-stats
    [hh,T] = evalc('topoplot_stats(cfg, summary.channels, result.tstats);');
    if params.general.show_warnings
        fprintf(T);
    end
    title(band);
    
end

saveas(h, sprintf('%s/eeg_hilbert_bl_overtake_topo.png', outdir));

if show_plots
    h.Visible = 'on';
end


%    1.1.3 Topoplots - Negative vs. Positive - T-stat by channel

% Load layout
cfg = [];
cfg.layout = params.eeg.layout;
evalc('layout=ft_prepare_layout(cfg);');
cfg.layout = layout;
cfg.colorbar = 'yes';
cfg.marker = 'no';
cfg.zlim = [-5 5]; % params.eeg.hilbert.plot.topo.baseline_overtake.clim;
cfg.colorbar = 'yes';
cfg.colorbartext = 'T-stat';
cfg.colormap = '-RdBu'; %params.eeg.hilbert.plot.topo.baseline_overtake.colormap;

h = figure('Visible','off');
h.Color = 'w';

h.Position([3 4]) = [320*N_bands 250];
t = tiledlayout("horizontal");
title(t,'Hilbert transform (Negative-Positive Outcome)');

for bb = 1 : N_bands
    nexttile;

    % Prepare stats
    band = summary.bands.Band{bb};
    result = summary.stats.ttest.channels.(band).overtake.outcomes{1};
    
    % Plot t-stats
    [hh,T] = evalc('topoplot_stats(cfg, summary.channels, result.tstats);');
    if params.general.show_warnings
        fprintf(T);
    end
    title(band);
    
end

saveas(h, sprintf('%s/eeg_hilbert_overtake_oputcome_topo.png', outdir));

if show_plots
    h.Visible = 'on';
end



% 1.2.1. Boxplots - Overtake Outcome (Negative vs. Positive)

h = plotlyfig('Visible','off');
h.data = cell(3,1);
Y = zeros(N_bands*N_subj,2);
X = cell(N_bands*N_subj,1);
G = [];
idx = 0;

% Boxplots groups by frequency band
for bb = 1 : N_bands
    band = summary.bands.Band{bb};
    result = summary.stats.ttest.grand_mean.(band).overtake.outcomes{1};
    
    Y(idx+1:idx+N_subj,:) = result.data;
    X(idx+1:idx+N_subj) = {band};

    G = result.labels;
    idx = idx + N_subj;
end

for i = 1 : 2
    h.data(i) = {struct( ...
                    'y', Y(:,i), ...
                    'x', {X}, ...
                    'name', G(i), ...
                    'type', 'box' ...
                    )};
end

% Significance asterisks
Y = ones(N_bands,1)*0.25;
X = {};
T = {};

for bb = 1 : N_bands
    band = summary.bands.Band{bb};
    result = summary.stats.ttest.grand_mean.(band).overtake.outcomes{1};
    if result.sig_fdr
        X = [X {band}];
        T = [T {'-*-'}];
    end
end

h.data(3) = {struct( ...
                'y', Y, ...
                'x', {X}, ...
                'name', 'p<0.05', ...
                'type', 'scatter', ...
                'mode', 'text', ...
                'textposition', 'top', ...
                'text', {T}, ...
                'showlegend', false, ...
                'textfont', struct(...
                      'family', 'Courier New, monospace', ...
                      'size', 40, ...
                      'color', 'black')...
                )};


h.layout = plotly_layouts.layouts.hilbert.boxplot;
h.layout.title = 'Hilbert envelopes: Negative v. Positive Outcome';

h.PlotOptions.SaveFolder = outdir;
h.PlotOptions.FileName = 'eeg_hilbert_overtake_outcome_boxplots';

try
    plotlyoffline(h);
catch ex
    fprintf('\nAttempting to fetch plotlyoffline dependencies...');
    getplotlyoffline(params.general.plotly_url);
    plotlyoffline(h);
end

if show_plots
    web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
end


%    1.3.1. Line plots - across cycles/rounds, average over channels

N_cycles = max(summary.N_cycles);
h = figure('Visible','off');
h.Color = 'w';
clr = lines(N_cycles);
offset = 0.5;

% Separate plots for each frequency band
for bb = 1 : N_bands

    band = summary.bands.Band{bb};
    M = summary.stats.glm.grand_mean.cycles.(band).means + offset*bb;
    S = summary.stats.glm.grand_mean.cycles.(band).stds;

    errorbar(1:N_cycles, M, S, 'Color', clr(bb,:), 'LineWidth', 2);

    hold on;

end

axis = gca;
axis.YAxis.TickValues = offset:offset:N_bands*offset;
axis.YAxis.TickLabels = summary.bands.Band;
axis.XAxis.TickValues = 1:N_cycles;
xlim([0 N_cycles+1]);
ylim([0 (N_bands+1)*offset]);

xlabel('Round');
axis.FontSize = 16;
title('Hilbert envelopes: Over Rounds');

filename = 'eeg_hilbert_rounds_lines';
fig2plotly(h, 'savefolder', outdir, ...
              'filename', 'eeg_hilbert_rounds_lines', ...
              'visible', false, ...
              'offline', true, ...
              'open', false);

if show_plots
    web(sprintf('file://%s/%s.html', outdir, filename), '-new', '-notoolbar');
end



% 1.4.1. Hilbert relationship to pupil diameter



end

