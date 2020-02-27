function [  ] = plot_event_summary_eye ( params, summary, out2file )
% Create plots from summary data
% This version uses the Plotly library
%
% Apparently this library cannot export as PNG yet...
%

if nargin < 3
   out2file = false;
end

show_plots = params.eye.epochs.plots.show_webplots;
plotly_layouts = load(params.general.plots.plotly_layouts_file);
sig_clrs = ['#7f7f7f' dec2hex(120)];
sig_dims = [-.25 .1];

outdir = sprintf('%s/figures', params.io.results_dir);
if ~exist(outdir, 'dir')
   mkdir(outdir); 
end

% % Difference from zero - Left change (Passing onset)
h = plot_timelocked(params, summary.stats.left_change, {'Passing Onset'}, ...
                    plotly_layouts.layouts.boxplot, sig_dims, sig_clrs);
h.layout.title = 'PD Locked to Passing Onset';
h.layout.yaxis.range = [-0.5 1.5]; % [min(y(:)) max(y(:))]*1.3;
h.PlotOptions.FileName = sprintf('%s/events_passonset_bl', outdir);
plotlyoffline(h);

if show_plots
    web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
end

% if out2file
%     saveplotlyfig(h, sprintf('%s/timelocked_passonset_bl.png', outdir));
% end

% Difference from zero - Right change (Passing offset)
h = plot_timelocked(params, summary.stats.right_change, {'Passing Offset'}, ...
                    plotly_layouts.layouts.boxplot, sig_dims, sig_clrs);
h.layout.title = 'PD Locked to Passing Offset';
h.layout.yaxis.range = [-0.5 1.5]; % [min(y(:)) max(y(:))]*1.3;
h.PlotOptions.FileName = sprintf('%s/events_passoffset_bl', outdir);
plotlyoffline(h);

if show_plots
    web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
end

% if out2file
%     saveplotlyfig(h, sprintf('%s/timelocked_passoffset_bl.png', outdir));
% end

% Difference from zero - Overtake event
h = plot_timelocked(params, summary.stats.overtake, {'Overtake'}, ...
                    plotly_layouts.layouts.boxplot, sig_dims, sig_clrs);
h.layout.title = 'PD Locked to Overtake Event';
h.layout.yaxis.range = [-0.5 1.5]; % [min(y(:)) max(y(:))]*1.3;
h.PlotOptions.FileName = sprintf('%s/events_overtake_bl', outdir);
plotlyoffline(h);

if show_plots
    web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
end

% if out2file
%     saveplotlyfig(h, sprintf('%s/timelocked_overtake_bl.png', outdir));
% end

%Comparison - Easy vs. Difficult - Passing Onset 

h = plot_timelocked(params, summary.stats.left_change.diff, [{'Easy'},{'Difficult'}], ...
                    plotly_layouts.layouts.boxplot, sig_dims, sig_clrs);
h.layout.title = 'PD Locked to Passing Onset: Difficulty';
h.layout.yaxis.range = [-0.5 1.5]; % [min(y(:)) max(y(:))]*1.3;
h.layout.showlegend = true;
h.layout.legend = struct('x', 0.1, 'y', 0.95);

h.PlotOptions.FileName = sprintf('%s/events_onset_diff', outdir);
plotlyoffline(h);

if show_plots
    web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
end

% if out2file
%     saveplotlyfig(h, sprintf('%s/timelocked_onset_diff.png', outdir));
% end

% Comparison - Easy vs. Difficult - Passing Offset

h = plot_timelocked(params, summary.stats.right_change.diff, [{'Easy'},{'Difficult'}], ...
                    plotly_layouts.layouts.boxplot, sig_dims, sig_clrs);
h.layout.title = 'PD Locked to Passing Offset: Difficulty';
h.layout.yaxis.range = [-0.5 1.5]; % [min(y(:)) max(y(:))]*1.3;
h.layout.showlegend = true;
h.layout.legend = struct('x', 0.1, 'y', 0.95);

h.PlotOptions.FileName = sprintf('%s/events_offset_diff', outdir);
plotlyoffline(h);

if show_plots
    web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
end

% if out2file
%     saveplotlyfig(h, sprintf('%s/timelocked_offset_diff.png', outdir));
% end

% Comparison - Easy vs. Difficult - Overtake Event

h = plot_timelocked(params, summary.stats.overtake.diff, [{'Easy'},{'Difficult'}], ...
                    plotly_layouts.layouts.boxplot, sig_dims, sig_clrs);
h.layout.title = 'PD Locked to Overtake Event: Difficulty';
h.layout.yaxis.range = [-0.5 1.5]; % [min(y(:)) max(y(:))]*1.3;
h.layout.showlegend = true;
h.layout.legend = struct('x', 0.1, 'y', 0.95);

h.PlotOptions.FileName = sprintf('%s/events_overtake_diff', outdir);
plotlyoffline(h);

if show_plots
    web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
end

% if out2file
%     saveplotlyfig(h, sprintf('%s/timelocked_overtake_diff.png', outdir));
% end


% Comparison - Positive vs. Negative Outcome - Passing Onset

h = plot_timelocked(params, summary.stats.left_change.outcomes, [{'Negative'},{'Positive'}], ...
                    plotly_layouts.layouts.boxplot, sig_dims, sig_clrs);
h.layout.title = 'PD Locked to Passing Onset: Outcome';
h.layout.yaxis.range = [-0.5 2]; % [min(y(:)) max(y(:))]*1.3;
h.layout.showlegend = true;
h.layout.legend = struct('x', 0.1, 'y', 0.95);

h.PlotOptions.FileName = sprintf('%s/events_onset_outcome', outdir);
plotlyoffline(h);

if show_plots
    web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
end

% if out2file
%     saveplotlyfig(h, sprintf('%s/timelocked_onset_outcome.png', outdir));
% end

% Comparison - Positive vs. Negative Outcome - Passing Offset

h = plot_timelocked(params, summary.stats.right_change.outcomes, [{'Negative'},{'Positive'}], ...
                    plotly_layouts.layouts.boxplot, sig_dims, sig_clrs);
h.layout.title = 'PD Locked to Passing Offset: Outcome';
h.layout.yaxis.range = [-0.5 2]; % [min(y(:)) max(y(:))]*1.3;
h.layout.showlegend = true;
h.layout.legend = struct('x', 0.1, 'y', 0.95);

h.PlotOptions.FileName = sprintf('%s/events_offset_outcome', outdir);
plotlyoffline(h);

if show_plots
    web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
end

% if out2file
%     saveplotlyfig(h, sprintf('%s/timelocked_offset_outcome.png', outdir));
% end

% Comparison - Positive vs. Negative Outcome - Overtake Event

h = plot_timelocked(params, summary.stats.overtake.outcomes, [{'Negative'},{'Positive'}], ...
                    plotly_layouts.layouts.boxplot, sig_dims, sig_clrs);
h.layout.title = 'PD Locked to Overtake Event: Outcome';
h.layout.yaxis.range = [-0.5 2]; % [min(y(:)) max(y(:))]*1.3;
h.layout.showlegend = true;
h.layout.legend = struct('x', 0.1, 'y', 0.95);

h.PlotOptions.FileName = sprintf('%s/events_overtake_outcome', outdir);
plotlyoffline(h);

if show_plots
    web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
end

% if out2file
%     saveplotlyfig(h, sprintf('%s/timelocked_overtake_outcome.png', outdir));
% end



end