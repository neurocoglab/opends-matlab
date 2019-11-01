function [ h ] = plot_epoch_summary2 ( summary, params, out2file, showplots )
% Create plots from summary data
% This version uses the Plotly library
%

%     '#1f77b4',  // muted blue
%     '#ff7f0e',  // safety orange
%     '#2ca02c',  // cooked asparagus green
%     '#d62728',  // brick red
%     '#9467bd',  // muted purple
%     '#8c564b',  // chestnut brown
%     '#e377c2',  // raspberry yogurt pink
%     '#7f7f7f',  // middle gray
%     '#bcbd22',  // curry yellow-green
%     '#17becf'   // blue-teal

plotly_clrs = [{'#1f77b4'},{'#ff7f0e'},{'#2ca02c'}];

if nargin < 4
   out2file = false;
end

if nargin < 5
    showplots = true;
end

outdir = sprintf('%s/plotly', params.output_dir);
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
h.layout = params.epochs.plots.layouts.boxplot;
h.layout.title = 'Baseline v. Passing Onset';
h.layout.yaxis.range = [min(min(y_baseline,y_passing)) max(max(y_baseline,y_passing))]*1.3;

h.PlotOptions.FileName = sprintf('%s/boxplots_pd_zscore_pass_baseline', outdir);
plotlyoffline(h);

if showplots
    web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
end

if out2file
    saveplotlyfig(h, sprintf('%s/boxplots_pd_zscore_pass_baseline.png', outdir));
end

% By difficulty
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
h.layout = params.epochs.plots.layouts.boxplot;
h.layout.title = 'Passing Onset: Easy v. Difficult';
h.layout.yaxis.range = [min(min(y_easy,y_diff)) max(max(y_easy,y_diff))]*1.3;

if out2file
    h.PlotOptions.FileName = sprintf('%s/boxplots_pd_zscore_pass_difficulty', outdir);
end
plotlyoffline(h);

if showplots
    web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
end

if out2file
    saveplotlyfig(h, sprintf('%s/boxplots_pd_zscore_pass_difficulty.png', outdir));
end


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
h.layout = params.epochs.plots.layouts.boxplot;
h.layout.title = 'Passing Onset: Positive v. Negative Outcome';
h.layout.yaxis.range = [min(min(y_pos,y_neg)) max(max(y_pos,y_neg))]*1.3;

if out2file
    h.PlotOptions.FileName = sprintf('%s/boxplots_pd_zscore_pass_outcome', outdir);
end
plotlyoffline(h);

if showplots
    web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
end

if out2file
    saveplotlyfig(h, sprintf('%s/boxplots_pd_zscore_pass_outcome.png', outdir));
end


if out2file || params.epochs.plot_scatter

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
    h.layout = params.epochs.plots.layouts.lines;
    h.layout.title = 'Baseline v. Passing By Subject';
    h.layout.yaxis.range = [min(min(y_baseline,y_passing)) max(max(y_baseline,y_passing))]*1.3;
    h.layout.xaxis.ticktext = [{'Baseline'},{'Passing'}];
    
    h.PlotOptions.FileName = sprintf('%s/lines_baseline_pass_subjects', outdir);
    plotlyoffline(h);
    
    if showplots
        web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
    end
    
    if out2file
        saveplotlyfig(h, sprintf('%s/lines_baseline_pass_subjects.png', outdir));
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
    yrange = [min(min(y_baseline,y_passing)) max(max(y_baseline,y_passing))]*1.3;
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
                    'name', 'Passing Onset', ...
                    'visible', 1, ...
                    'mode', 'lines', ...
                    'line', struct('color', plotly_clrs{2}, ...
                                   'width', 5) ...
                    )};
    
    h = plotlyfig('Visible','off');
    h.data = [data data2];
    h.layout = params.epochs.plots.layouts.boxplot;
    h.layout.title = 'PD v. Simulation Score';
    h.layout.yaxis.range = yrange;
    h.layout.xaxis.title = 'Simulation Score';
    h.layout.showlegend = true;
    h.layout.legend = struct('x', 0.05, 'y', 0.05);
    
    h.PlotOptions.FileName = sprintf('%s/scatter_pdz_X_score', outdir);
    plotlyoffline(h);
    
    if showplots
        web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
    end
    
    if out2file
        saveplotlyfig(h, sprintf('%s/scatter_pdz_X_score.png', outdir));
    end
    
    % Plot Difference PD (Overtake - Baseline) versus Score
    if out2file
        h = figure('visible','off');
    else
        h = figure;
    end 
    
    lm1 = fitlm(scores,diff_bp);
    y1 = [lm1.predict(xrange(1)) lm1.predict(xrange(2))];
    
    data = {struct(...
                    'x', scores, ...
                    'y', diff_bp, ...
                    'name', 'PassingOnset-Baseline', ...
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
    
   
    yrange = [min(diff_bp) max(diff_bp)]*1.3;
    h = plotlyfig('Visible','off');
    h.data = [data];
    h.layout = params.epochs.plots.layouts.boxplot;
    h.layout.title = 'PD (Passing Onset minus Baseline) v. Simulation Score';
    h.layout.yaxis.range = yrange;
    h.layout.xaxis.title = 'Simulation Score';
    h.layout.yaxis.title = 'Mean PD Difference (Z-score)';
    
    h.PlotOptions.FileName = sprintf('%s/scatter_diff_pdz_X_score', outdir);
    plotlyoffline(h);
    
    if showplots
        web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
    end
    
    if out2file
        saveplotlyfig(h, sprintf('%s/scatter_diff_pdz_X_score.png', outdir));
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

                                    
if out2file
    h = figure('visible','off');
else
    h = figure;
end

scores = [1:N_cycles; 1:N_cycles]';
y = [pd_bl_mean,pd_pass_mean];
lower = [pd_bl_ci(:,1),pd_pass_ci(:,1)]; % / N_subs;
upper = [pd_bl_ci(:,2),pd_pass_ci(:,2)];

h = plotlyfig('Visible','off');
[data_line,data_shade] = get_plotly_ci_data(scores, y, lower, upper, ...
                                            plotly_clrs, [{'Baseline'},{'Passing Onset'}], ...
                                            0.1);
h.data = [data_line;data_shade];
h.layout = params.epochs.plots.layouts.boxplot;
h.layout.title = 'PD Across Simulation Rounds';
h.layout.yaxis.range = [min(y(:)) max(y(:))]*1.3;
h.layout.xaxis.title = 'Simulation Round';
h.layout.showlegend = true;
h.layout.legend = struct('x', 0.75, 'y', 0.95);

h.PlotOptions.FileName = sprintf('%s/lines_pdz_rounds', outdir);
plotlyoffline(h);

if showplots
    web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
end

if out2file
    saveplotlyfig(h, sprintf('%s/lines_pdz_rounds.png', outdir));
end



% Outcomes v. Difficulty

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
h.layout = params.epochs.plots.layouts.boxplot;
h.layout.boxmode = 'group';
h.layout.yaxis.range = [min(y) max(y)] * 1.3;
h.layout.yaxis.zeroline = false;
h.layout.boxgroupgap = 0.5;
h.layout.boxgap = 0.2;

h.layout.title = 'PD by Difficulty and Outcome';
% h.layout.xaxis.title = 'Difficulty';
h.layout.showlegend = true;

h.PlotOptions.FileName = sprintf('%s/boxplots_pdz_diff_outcome', outdir);
plotlyoffline(h);

if showplots
    web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
end

if out2file
    saveplotlyfig(h, sprintf('%s/boxplots_pdz_diff_outcome.png', outdir));
end

end