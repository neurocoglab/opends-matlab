function [ h ] = plot_epoch_summary_eye ( params, summary, out2file )
% Create plots from summary data
% This version uses the Plotly library
%

if nargin < 3
   out2file = false;
end

show_plots = params.eye.epochs.plots.show_webplots;
plotly_clrs = params.general.plots.plotly_colors;
plotly_layouts = load(params.general.plots.plotly_layouts_file);

outdir = sprintf('%s/figures', params.io.results_dir);
if ~exist(outdir, 'dir')
   mkdir(outdir); 
end

% Full session (z-score)
tbl = summary.stats.passing_baseline.data;

y_baseline = tbl.PD(tbl.IsBaseline==1,:);
y_passing = tbl.PD(tbl.IsBaseline==0,:);

data = {...
  struct(...
    'y', y_baseline, ...
    'name', 'Baseline', ...
    'boxpoints', 'all', ...
    'jitter', 0.3, ...
    'pointpos', -1.8, ...
    'type', 'box'), ...
  struct(...
    'y', y_passing, ...
    'name', 'Passing', ...
    'boxpoints', 'all', ...
    'jitter', 0.3, ...
    'pointpos', -1.8, ...
    'type', 'box')
};

h = plotlyfig('Visible','off'); % initalize an empty figure object
h.data = data;
h.layout = plotly_layouts.layouts.boxplot;
h.layout.title = 'Baseline v. Passing Epochs';
yrange = [min(min(y_baseline,y_passing)) max(max(y_baseline,y_passing))];
ypad = abs(diff(yrange))*0.3/2;
yrange(1) = yrange(1)-ypad;yrange(2) = yrange(2)+ypad;
h.layout.yaxis.range = yrange;

h.PlotOptions.SaveFolder = outdir;
h.PlotOptions.FileName = 'epochs_boxplots_pd_zscore_pass_baseline'; % sprintf('%s/epochs_boxplots_pd_zscore_pass_baseline', outdir);
plotlyoffline(h);

if show_plots
    web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
end

% if out2file
%     saveplotlyfig(h, sprintf('%s/boxplots_pd_zscore_pass_baseline.png', outdir));
% end

% By difficulty
if params.sim.epochs.difficulty.apply
    tbl = summary.stats.passing_outcome_diff.data;

    y_easy = tbl.PD(strcmp(tbl.Difficulty,'Easy'),:);
    y_diff = tbl.PD(strcmp(tbl.Difficulty,'Difficult'),:);

    data = {...
      struct(...
        'y', y_easy, ...
        'name', 'Easy', ...
        'boxpoints', 'all', ...
        'jitter', 0.3, ...
        'pointpos', -1.8, ...
        'type', 'box'), ...
      struct(...
        'y', y_diff, ...
        'name', 'Difficult', ...
        'boxpoints', 'all', ...
        'jitter', 0.3, ...
        'pointpos', -1.8, ...
        'type', 'box')
    };

    h = plotlyfig('Visible','off'); % initalize an empty figure object
    h.data = data;
    h.layout = plotly_layouts.layouts.boxplot;
    h.layout.title = 'Passing Epochs: Easy v. Difficult';
    yrange = [min(min(y_easy,y_diff)) max(max(y_easy,y_diff))];
    ypad = abs(diff(yrange))*0.3/2;
    yrange(1) = yrange(1)-ypad;yrange(2) = yrange(2)+ypad;
    h.layout.yaxis.range = yrange;

    h.PlotOptions.SaveFolder = outdir;
    h.PlotOptions.FileName = 'epochs_boxplots_pd_zscore_pass_difficulty';
    plotlyoffline(h);

    if show_plots
        web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
    end
else
    tbl = summary.stats.passing_outcome.data;
end

% if out2file
%     saveplotlyfig(h, sprintf('%s/epochs_boxplots_pd_zscore_pass_difficulty.png', outdir));
% end


% By outcome

y_pos = tbl.PD(strcmp(tbl.Outcome,'Positive'),:);
y_neg = tbl.PD(strcmp(tbl.Outcome,'Negative'),:);

data = {...
  struct(...
    'y', y_pos, ...
    'name', 'Positive', ...
    'boxpoints', 'all', ...
    'jitter', 0.3, ...
    'pointpos', -1.8, ...
    'type', 'box'), ...
  struct(...
    'y', y_neg, ...
    'name', 'Negative', ...
    'boxpoints', 'all', ...
    'jitter', 0.3, ...
    'pointpos', -1.8, ...
    'type', 'box')
};

h = plotlyfig('Visible','off'); % initalize an empty figure object
h.data = data;
h.layout = plotly_layouts.layouts.boxplot;
h.layout.title = 'Passing Epochs: Positive v. Negative Outcome';
yrange = [min(min(y_pos,y_neg)) max(max(y_pos,y_neg))];
ypad = abs(diff(yrange))*0.3/2;
yrange(1) = yrange(1)-ypad;yrange(2) = yrange(2)+ypad;
h.layout.yaxis.range = yrange;

h.PlotOptions.SaveFolder = outdir;
h.PlotOptions.FileName = 'epochs_boxplots_pd_zscore_pass_outcome';

plotlyoffline(h);

if show_plots
    web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
end

% if out2file
%     saveplotlyfig(h, sprintf('%s/boxplots_pd_zscore_pass_outcome.png', outdir));
% end


if out2file || params.eye.epochs.scatter

    % Subject score versus pupil dilation
    subjects = summary.stats.passing_baseline_dscore.data.Subject;
    scores = summary.stats.passing_baseline_dscore.data.Scores;
    diff_bp = summary.stats.passing_baseline_dscore.data.DeltaPD;
    y_baseline = summary.stats.passing_baseline_dscore.data.PD_baseline;
    y_passing = summary.stats.passing_baseline_dscore.data.PD_passing;

    clr = [0.4 1 1 1];
    mi=min(scores)-200; mx = max(scores);
    data = {};

    for i = 1 : length(subjects)
        cc = 1 - clr * (scores(i)-mi)/(mx-mi);
        cc = round(cc * 255);
        data(end+1) = {struct(...
                        'x', [1 2], ...
                        'y', [y_baseline(i) y_passing(i)], ...
                        'color', 'rgb(200,0,0)', ...
                        'name', subjects{i}, ...
                        'visible', 1, ...
                        'mode', 'lines+markers', ...
                        'marker', struct(...
                                        'symbol', 'circle', ...
                                        'size', 13, ...
                                        'color', sprintf('rgba(%d,%d,%d,%d)', cc(1), cc(2), cc(3), cc(4)) ...
                                        ) ...
                              ) ...
                      };
    end

    h = plotlyfig('Visible','off'); % initalize an empty figure object
    h.data = data;
    h.layout = plotly_layouts.layouts.lines;
    h.layout.title = 'Baseline v. Passing By Subject';
    
    yrange = [min(min(y_baseline,y_passing)) max(max(y_baseline,y_passing))];
    ypad = abs(diff(yrange))*0.3/2;
    yrange(1) = yrange(1)-ypad;yrange(2) = yrange(2)+ypad;
    h.layout.yaxis.range = yrange;
    h.layout.xaxis.ticktext = [{'Baseline'},{'Passing'}];
    
    h.PlotOptions.SaveFolder = outdir;
    h.PlotOptions.FileName = 'epochs_lines_baseline_pass_subjects';
    plotlyoffline(h);
    
    if show_plots
        web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
    end
    
    % Scatterplots (z-score)
    % Plot Overtake & Baseline PD versus Score

    data = {struct(...
                    'x', scores, ...
                    'y', y_baseline, ...
                    'name', 'Baseline', ...
                    'visible', 1, ...
                    'mode', 'markers', ...
                    'marker', struct(...
                                    'symbol', 'circle', ...
                                    'size', 13 ...
                                    ) ...
                    ), ...
            struct(...
                    'x', scores, ...
                    'y', y_passing, ...
                    'name', 'Passing', ...
                    'visible', 1, ...
                    'mode', 'markers', ...
                    'marker', struct(...
                                    'symbol', 'circle', ...
                                    'size', 13 ...
                                    ) ...
                    )};
                
    % Add regression lines
    yrange = [min(min(y_baseline,y_passing)) max(max(y_baseline,y_passing))];
    ypad = abs(diff(yrange))*0.3/2;
    yrange(1) = yrange(1)-ypad;yrange(2) = yrange(2)+ypad;
    h.layout.yaxis.range = yrange;
    xrange = [min(scores) max(scores)];
    lm1 = fitlm(scores,y_baseline);
    lm2 = fitlm(scores,y_passing);
    y1 = [lm1.predict(xrange(1)) lm1.predict(xrange(2))];
    y2 = [lm2.predict(xrange(1)) lm2.predict(xrange(2))];
    data2 = {struct(...
                    'showlegend', false, ...
                    'x', xrange, ...
                    'y', y1, ...
                    'name', 'Baseline', ...
                    'visible', 1, ...
                    'mode', 'lines', ...
                    'line', struct('color', plotly_clrs{1}, ...
                                   'width', 5) ...
                    ), ...
            struct(...
                    'showlegend', false, ...
                    'x', xrange, ...
                    'y', y2, ...
                    'name', 'Passing', ...
                    'visible', 1, ...
                    'mode', 'lines', ...
                    'line', struct('color', plotly_clrs{2}, ...
                                   'width', 5) ...
                    )};
    
    h = plotlyfig('Visible','off');
    h.data = [data data2];
    h.layout = plotly_layouts.layouts.boxplot;
    h.layout.title = 'PD v. Simulation Score';
    h.layout.yaxis.range = yrange;
    h.layout.xaxis.title = 'Simulation Score';
    h.layout.showlegend = true;
    h.layout.legend = struct('x', 0.05, 'y', 0.05);
    
    h.PlotOptions.SaveFolder = outdir;
    h.PlotOptions.FileName = 'epochs_scatter_pdz_X_score';
    plotlyoffline(h);
    
    if show_plots
        web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
    end
    
    
    lm1 = fitlm(scores,diff_bp);
    y1 = [lm1.predict(xrange(1)) lm1.predict(xrange(2))];
    
    data = {struct(...
                    'x', scores, ...
                    'y', diff_bp, ...
                    'name', 'Passing-Baseline', ...
                    'visible', 1, ...
                    'mode', 'markers', ...
                    'marker', struct(...
                                    'symbol', 'circle', ...
                                    'size', 13 ...
                                    ) ...
                    ), ...
             struct(...
                    'x', xrange, ...
                    'y', y1, ...
                    'name', 'Regression Line', ...
                    'visible', 1, ...
                    'mode', 'lines', ...
                    'line', struct('color', plotly_clrs{1}, ...
                                   'width', 5) ...
                    )};
    
    yrange = [min(diff_bp) max(diff_bp)];
    ypad = abs(diff(yrange))*0.3/2;
    yrange(1) = yrange(1)-ypad;yrange(2) = yrange(2)+ypad;
    h.layout.yaxis.range = yrange;
    h = plotlyfig('Visible','off');
    h.data = [data];
    h.layout = plotly_layouts.layouts.boxplot;
    h.layout.title = 'PD (Passing minus Baseline Epochs) v. Simulation Score';
    h.layout.yaxis.range = yrange;
    h.layout.xaxis.title = 'Simulation Score';
    h.layout.yaxis.title = 'Mean PD Difference (Z-score)';
    
    h.PlotOptions.SaveFolder = outdir;
    h.PlotOptions.FileName = 'epochs_scatter_diff_pdz_X_score';
    plotlyoffline(h);
    
    if show_plots
        web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
    end

end

% Plot cycles
N_cycles = length(summary.cycles.baseline.pupil);

% Baseline - Boxplots


% Passing - Boxplots


% Baseline & Passing - Line plot

pd_bl_mean = zeros(N_cycles,1);
pd_bl_ci = zeros(N_cycles,2);
pd_pass_mean = zeros(N_cycles,1);
pd_pass_ci = zeros(N_cycles,2);

N_subs = length(summary.subjects);

pd_bl = zeros(N_cycles, N_subs);
pd_pass = zeros(N_cycles, N_subs);

for j = 1 : N_cycles
    x_bl = zeros(N_subs,1);
    x_pass = zeros(N_subs,1);
    for i = 1 : N_subs
       if length(summary.zscore.cycles.baseline.subjects.pupil{i}) < N_cycles
           %ffs
           x_bl(i) = nan;
           x_pass(i) = nan;
       else
           x_bl(i) = mean(summary.zscore.cycles.baseline.subjects.pupil{i}{j});
           x_pass(i) = mean(summary.zscore.cycles.passing.subjects.pupil{i}{j});
           
       end
    end
    
    pd_bl(j,:) = x_bl;
    pd_pass(j,:) = x_pass;
    
    pd_bl_mean(j) = nanmean(x_bl);
    pd_pass_mean(j) = nanmean(x_pass);
    
    pdist = fitdist(x_bl(~isnan(x_bl)),'Normal');
    ci = paramci(pdist);
    pd_bl_ci(j,:) = ci(:,1);
    pdist = fitdist(x_pass(~isnan(x_pass)),'Normal');
    ci = paramci(pdist);
    pd_pass_ci(j,:) = ci(:,1);
    
end

scores = [1:N_cycles; 1:N_cycles]';
y = [pd_bl_mean,pd_pass_mean];
lower = [pd_bl_ci(:,1),pd_pass_ci(:,1)]; % / N_subs;
upper = [pd_bl_ci(:,2),pd_pass_ci(:,2)];

h = plotlyfig('Visible','off');
[data_line,data_shade] = get_plotly_ci_data(scores, y, lower, upper, ...
                                            [{'Baseline'},{'Passing'}], plotly_clrs, ...
                                            0.1);
h.data = [data_line;data_shade];
h.layout = plotly_layouts.layouts.boxplot;
h.layout.title = 'PD Across Simulation Rounds';
yrange = [min(lower(:)) max(upper(:))];
ypad = abs(diff(yrange))*0.3/2;
yrange(1) = yrange(1)-ypad;yrange(2) = yrange(2)+ypad;
h.layout.yaxis.range = yrange;
h.layout.xaxis.title = 'Simulation Round';
h.layout.showlegend = true;
h.layout.legend = struct('x', 0.75, 'y', 0.95);

h.PlotOptions.SaveFolder = outdir;
h.PlotOptions.FileName = 'epochs_lines_pdz_rounds';
plotlyoffline(h);

if show_plots
    web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
end

% Score by round



% Outcomes v. Difficulty
if params.sim.epochs.difficulty.apply
    tbl = summary.stats.passing_outcome_diff.data;
    idx_diff_pos = find(strcmp(tbl.Difficulty,'Difficult') & strcmp(tbl.Outcome,'Positive'));
    idx_diff_neg = find(strcmp(tbl.Difficulty,'Difficult') & strcmp(tbl.Outcome,'Negative'));
    idx_easy_pos = find(strcmp(tbl.Difficulty,'Easy') & strcmp(tbl.Outcome,'Positive'));
    idx_easy_neg = find(strcmp(tbl.Difficulty,'Easy') & strcmp(tbl.Outcome,'Negative'));

    x = [repmat({'Difficult'}, 1, length(idx_diff_pos)),repmat({'Easy'}, 1, length(idx_easy_pos))];
    % x = [ones(length(idx_diff_pos),1);zeros(length(idx_easy_pos),1)]';
    y = [tbl.PD(idx_diff_pos);tbl.PD(idx_easy_pos)];
    data1 = struct(...
        'x', {x'}, ...
        'y', y, ...
        'name', 'Positive', ...
        'boxpoints', 'all', ...
        'jitter', 0.3, ...
        'pointpos', -1.8, ...
        'type', 'box');

    maxy=max(y);

    x = [repmat({'Difficult'}, 1, length(idx_diff_neg)),repmat({'Easy'}, 1, length(idx_easy_neg))];
    % x = [ones(length(idx_diff_pos),1);zeros(length(idx_easy_pos),1)]';
    y = [tbl.PD(idx_diff_neg);tbl.PD(idx_easy_neg)];
    data2 = struct(...
        'x', {x'}, ...
        'y', y, ...
        'name', 'Negative', ...
        'boxpoints', 'all', ...
        'jitter', 0.3, ...
        'pointpos', -1.8, ...
        'type', 'box');

    h = plotlyfig('Visible','off');

    h.data = {data1,data2};
    h.layout = plotly_layouts.layouts.boxplot;
    h.layout.boxmode = 'group';
    yrange = [min(y) max(maxy, max(y))];
    ypad = abs(diff(yrange))*0.3/2;
    yrange(1) = yrange(1)-ypad;yrange(2) = yrange(2)+ypad;
    h.layout.yaxis.range = yrange;
    h.layout.yaxis.zeroline = false;
    h.layout.boxgroupgap = 0.5;
    h.layout.boxgap = 0.2;

    h.layout.title = 'PD by Difficulty and Outcome';
    % h.layout.xaxis.title = 'Difficulty';
    h.layout.showlegend = true;

    h.PlotOptions.SaveFolder = outdir;
    h.PlotOptions.FileName = 'epochs_boxplots_pdz_diff_outcome';
    plotlyoffline(h);

    if show_plots
        web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
    end

end



end