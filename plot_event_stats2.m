function [] = plot_event_stats2 ( summary, params, out2file, showplots )
% Create plots from summary data
% This version uses the Plotly library
%

if nargin < 3
   out2file = false;
end

if nargin < 4
    showplots = true;
end

outdir = sprintf('%s/plotly', params.output_dir);
if ~exist(outdir, 'dir')
   mkdir(outdir); 
end

% % Difference from zero - Left change (Passing onset)
N_subj = length(summary.subjects);

h = plot_timelocked(summary.stats.left_change, N_subj, {'Passing Onset'}, params.epochs.plots.layouts.boxplot, [-.25 .1], ['#7f7f7f' dec2hex(120)]);
h.layout.title = 'PD Locked to Passing Onset';
h.layout.yaxis.range = [-0.5 1.5]; % [min(y(:)) max(y(:))]*1.3;
h.PlotOptions.FileName = sprintf('%s/timelocked_passonset_bl', outdir);
plotlyoffline(h);

if showplots
    web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
end

if out2file
    saveplotlyfig(h, sprintf('%s/timelocked_passonset_bl.png', outdir));
end

% Difference from zero - Right change (Passing offset)
h = plot_timelocked(summary.stats.right_change, N_subj, {'Passing Offset'}, params.epochs.plots.layouts.boxplot, [-.25 .1], ['#7f7f7f' dec2hex(120)]);
h.layout.title = 'PD Locked to Passing Offset';
h.layout.yaxis.range = [-0.5 1.5]; % [min(y(:)) max(y(:))]*1.3;
h.PlotOptions.FileName = sprintf('%s/timelocked_passoffset_bl', outdir);
plotlyoffline(h);

if showplots
    web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
end

if out2file
    saveplotlyfig(h, sprintf('%s/timelocked_passoffset_bl.png', outdir));
end

% Difference from zero - Overtake event
h = plot_timelocked(summary.stats.overtake, N_subj, {'Overtake'}, params.epochs.plots.layouts.boxplot, [-.25 .1], ['#7f7f7f' dec2hex(120)]);
h.layout.title = 'PD Locked to Overtake Event';
h.layout.yaxis.range = [-0.5 1.5]; % [min(y(:)) max(y(:))]*1.3;
h.PlotOptions.FileName = sprintf('%s/timelocked_overtake_bl', outdir);
plotlyoffline(h);

if showplots
    web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
end

if out2file
    saveplotlyfig(h, sprintf('%s/timelocked_overtake_bl.png', outdir));
end

%Comparison - Easy vs. Difficult - Passing Onset 

h = plot_timelocked(summary.stats.left_change.diff, N_subj, [{'Easy'},{'Difficult'}], params.epochs.plots.layouts.boxplot, [-.25 .1], ['#7f7f7f' dec2hex(120)]);
h.layout.title = 'PD Locked to Passing Onset: Difficulty';
h.layout.yaxis.range = [-0.5 1.5]; % [min(y(:)) max(y(:))]*1.3;
h.layout.showlegend = true;
h.layout.legend = struct('x', 0.1, 'y', 0.95);

h.PlotOptions.FileName = sprintf('%s/timelocked_onset_diff', outdir);
plotlyoffline(h);

if showplots
    web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
end

if out2file
    saveplotlyfig(h, sprintf('%s/timelocked_onset_diff.png', outdir));
end

% Comparison - Easy vs. Difficult - Passing Offset

h = plot_timelocked(summary.stats.right_change.diff, N_subj, [{'Easy'},{'Difficult'}], params.epochs.plots.layouts.boxplot, [-.25 .1], ['#7f7f7f' dec2hex(120)]);
h.layout.title = 'PD Locked to Passing Offset: Difficulty';
h.layout.yaxis.range = [-0.5 1.5]; % [min(y(:)) max(y(:))]*1.3;
h.layout.showlegend = true;
h.layout.legend = struct('x', 0.1, 'y', 0.95);

h.PlotOptions.FileName = sprintf('%s/timelocked_offset_diff', outdir);
plotlyoffline(h);

if showplots
    web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
end

if out2file
    saveplotlyfig(h, sprintf('%s/timelocked_offset_diff.png', outdir));
end

% Comparison - Easy vs. Difficult - Overtake Event

h = plot_timelocked(summary.stats.overtake.diff, N_subj, [{'Easy'},{'Difficult'}], params.epochs.plots.layouts.boxplot, [-.25 .1], ['#7f7f7f' dec2hex(120)]);
h.layout.title = 'PD Locked to Overtake Event: Difficulty';
h.layout.yaxis.range = [-0.5 1.5]; % [min(y(:)) max(y(:))]*1.3;
h.layout.showlegend = true;
h.layout.legend = struct('x', 0.1, 'y', 0.95);

h.PlotOptions.FileName = sprintf('%s/timelocked_overtake_diff', outdir);
plotlyoffline(h);

if showplots
    web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
end

if out2file
    saveplotlyfig(h, sprintf('%s/timelocked_overtake_diff.png', outdir));
end


% Comparison - Positive vs. Negative Outcome - Passing Onset

h = plot_timelocked(summary.stats.left_change.outcomes, N_subj, [{'Negative'},{'Positive'}], params.epochs.plots.layouts.boxplot, [-.25 .1], ['#7f7f7f' dec2hex(120)]);
h.layout.title = 'PD Locked to Passing Onset: Outcome';
h.layout.yaxis.range = [-0.5 2]; % [min(y(:)) max(y(:))]*1.3;
h.layout.showlegend = true;
h.layout.legend = struct('x', 0.1, 'y', 0.95);

h.PlotOptions.FileName = sprintf('%s/timelocked_onset_outcome', outdir);
plotlyoffline(h);

if showplots
    web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
end

if out2file
    saveplotlyfig(h, sprintf('%s/timelocked_onset_outcome.png', outdir));
end

% Comparison - Positive vs. Negative Outcome - Passing Offset

h = plot_timelocked(summary.stats.right_change.outcomes, N_subj, [{'Negative'},{'Positive'}], params.epochs.plots.layouts.boxplot, [-.25 .1], ['#7f7f7f' dec2hex(120)]);
h.layout.title = 'PD Locked to Passing Offset: Outcome';
h.layout.yaxis.range = [-0.5 2]; % [min(y(:)) max(y(:))]*1.3;
h.layout.showlegend = true;
h.layout.legend = struct('x', 0.1, 'y', 0.95);

h.PlotOptions.FileName = sprintf('%s/timelocked_offset_outcome', outdir);
plotlyoffline(h);

if showplots
    web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
end

if out2file
    saveplotlyfig(h, sprintf('%s/timelocked_offset_outcome.png', outdir));
end

% Comparison - Positive vs. Negative Outcome - Overtake Event

h = plot_timelocked(summary.stats.overtake.outcomes, N_subj, [{'Negative'},{'Positive'}], params.epochs.plots.layouts.boxplot, [-.25 .1], ['#7f7f7f' dec2hex(120)]);
h.layout.title = 'PD Locked to Overtake Event: Outcome';
h.layout.yaxis.range = [-0.5 2]; % [min(y(:)) max(y(:))]*1.3;
h.layout.showlegend = true;
h.layout.legend = struct('x', 0.1, 'y', 0.95);

h.PlotOptions.FileName = sprintf('%s/timelocked_overtake_outcome', outdir);
plotlyoffline(h);

if showplots
    web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
end

if out2file
    saveplotlyfig(h, sprintf('%s/timelocked_overtake_outcome.png', outdir));
end




    


end