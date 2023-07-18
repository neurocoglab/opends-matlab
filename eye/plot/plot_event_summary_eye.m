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
h.layout.yaxis.range = params.eye.events.right_change.plots.ylims; % [-0.5 1.5]; % [min(y(:)) max(y(:))]*1.3;
h.PlotOptions.SaveFolder = outdir;
h.PlotOptions.FileName = 'eye_events_passonset_bl';
plotlyoffline(h);

if show_plots
    web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
end

% if out2file
%     saveplotlyfig(h, sprintf('%s/timelocked_passonset_bl.png', outdir));
% end

% Difference from zero - Right change (Passing offset)
h = plot_timelocked(params, summary.stats.right_change, {'Passing Offset'}, ...
                    plotly_layouts.layouts.boxplot, sig_dims, sig_clrs);
h.layout.title = 'PD Locked to Passing Offset';
h.layout.yaxis.range = params.eye.events.right_change.plots.ylims; % [-0.5 1.5]; % [min(y(:)) max(y(:))]*1.3;
h.PlotOptions.SaveFolder = outdir;
h.PlotOptions.FileName = 'eye_events_passoffset_bl';
plotlyoffline(h);

if show_plots
    web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
end

% if out2file
%     saveplotlyfig(h, sprintf('%s/timelocked_passoffset_bl.png', outdir));
% end

% Difference from zero - Overtake event
h = plot_timelocked(params, summary.stats.overtake, {'Overtake'}, ...
                    plotly_layouts.layouts.boxplot, sig_dims, sig_clrs);
h.layout.title = 'PD Locked to Overtake Event';
h.layout.yaxis.range = params.eye.events.overtake.plots.ylims; % [-0.5 1.5]; % [min(y(:)) max(y(:))]*1.3;
h.PlotOptions.SaveFolder = outdir;
h.PlotOptions.FileName = 'eye_events_overtake_bl';
plotlyoffline(h);

if show_plots
    web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
end

% if out2file
%     saveplotlyfig(h, sprintf('%s/timelocked_overtake_bl.svg', outdir));
% end

% Comparison - Easy vs. Difficult 
if params.eye.events.difficulty.apply
    
    % Time series plot - Onset
    h = plot_timelocked(params, summary.stats.left_change.diff, [{'Easy'},{'Difficult'}], ...
                        plotly_layouts.layouts.boxplot, sig_dims, sig_clrs);
    h.layout.title = 'PD Locked to Passing Onset: Difficulty';
    h.layout.yaxis.range = params.eye.events.left_change.plots.ylims; [-0.5 1.5]; % [min(y(:)) max(y(:))]*1.3;
    h.layout.showlegend = true;
    h.layout.legend = struct('x', 0.1, 'y', 0.95);

    h.PlotOptions.SaveFolder = outdir;
    h.PlotOptions.FileName = 'eye_events_onset_diff';
    plotlyoffline(h);

    if show_plots
        web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
    end

    % Violin plot parameters (slope/amplitude) - Onset
    h.layout.title = 'PD Parameters for Passing Onset: Difficulty';
    
    h = plot_timelocked_params(params, summary.stats.left_change.diff.parameter_stats, ...
                              {'Easy','Difficult'});
    
    h.PlotOptions.SaveFolder = outdir;
    h.PlotOptions.FileName = 'eye_events_onset_diff_params';
    plotlyoffline(h);

    if show_plots
        web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
    end
    

    % Time series plot - Offset
    h = plot_timelocked(params, summary.stats.right_change.diff, [{'Easy'},{'Difficult'}], ...
                        plotly_layouts.layouts.boxplot, sig_dims, sig_clrs);
    h.layout.title = 'PD Locked to Passing Offset: Difficulty';
    h.layout.yaxis.range = params.eye.events.right_change.plots.ylims; % [-0.5 1.5]; % [min(y(:)) max(y(:))]*1.3;
    h.layout.showlegend = true;
    h.layout.legend = struct('x', 0.1, 'y', 0.95);

    h.PlotOptions.SaveFolder = outdir;
    h.PlotOptions.FileName = 'eye_events_offset_diff';
    plotlyoffline(h);

    if show_plots
        web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
    end

    % Time series plot - Overtake
    h = plot_timelocked(params, summary.stats.overtake.diff, [{'Easy'},{'Difficult'}], ...
                        plotly_layouts.layouts.boxplot, sig_dims, sig_clrs);
    h.layout.title = 'PD Locked to Overtake Event: Difficulty';
    h.layout.yaxis.range = params.eye.events.overtake.plots.ylims; % [-0.5 1.5]; % [min(y(:)) max(y(:))]*1.3;
    h.layout.showlegend = true;
    h.layout.legend = struct('x', 0.1, 'y', 0.95);

    h.PlotOptions.SaveFolder = outdir;
    h.PlotOptions.FileName = 'eye_events_overtake_diff';
    plotlyoffline(h);

    if show_plots
        web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
    end

end

% Comparison - Positive vs. Negative Outcome - Passing Onset
if params.eye.events.outcomes.apply
    h = plot_timelocked(params, summary.stats.left_change.outcomes, [{'Negative'},{'Positive'}], ...
                        plotly_layouts.layouts.boxplot, sig_dims, sig_clrs);
    h.layout.title = 'PD Locked to Passing Onset: Outcome';
    h.layout.yaxis.range = params.eye.events.left_change.plots.ylims; % [-0.5 2]; % [min(y(:)) max(y(:))]*1.3;
    h.layout.showlegend = true;
    h.layout.legend = struct('x', 0.1, 'y', 0.95);

     h.PlotOptions.SaveFolder = outdir;
    h.PlotOptions.FileName = 'eye_events_onset_outcome';
    plotlyoffline(h);

    if show_plots
        web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
    end

    % if out2file
    %     saveplotlyfig(h, sprintf('%s/timelocked_onset_outcome.png', outdir));
    % end

    % Comparison - Positive vs. Negative Outcome - Passing Offset

    h = plot_timelocked(params, summary.stats.right_change.outcomes, [{'Negative'},{'Positive'}], ...
                        plotly_layouts.layouts.boxplot, sig_dims, sig_clrs);
    h.layout.title = 'PD Locked to Passing Offset: Outcome';
    h.layout.yaxis.range = params.eye.events.right_change.plots.ylims; % [-0.5 2]; % [min(y(:)) max(y(:))]*1.3;
    h.layout.showlegend = true;
    h.layout.legend = struct('x', 0.1, 'y', 0.95);

    h.PlotOptions.SaveFolder = outdir;
    h.PlotOptions.FileName = 'eye_events_offset_outcome';
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
    h.layout.yaxis.range = params.eye.events.overtake.plots.ylims; % [-0.5 2]; % [min(y(:)) max(y(:))]*1.3;
    h.layout.showlegend = true;
    h.layout.legend = struct('x', 0.1, 'y', 0.95);

    h.PlotOptions.SaveFolder = outdir;
    h.PlotOptions.FileName = 'eye_events_overtake_outcome';
    plotlyoffline(h);

    if show_plots
        web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
    end

    % if out2file
    %     saveplotlyfig(h, sprintf('%s/timelocked_overtake_outcome.png', outdir));
    % end
end


%% Traffic events
if params.sim.events.traffic_decision.apply
    
    h = plot_timelocked(params, summary.stats.traffic_decision, {'Traffic Decision'}, ...
                        plotly_layouts.layouts.boxplot, sig_dims, sig_clrs);
    h.layout.title = 'PD Locked to Traffic Decision';
    h.layout.yaxis.range = params.eye.events.traffic_decision.plots.ylims; % [-0.5 2]; % [min(y(:)) max(y(:))]*1.3;
    h.layout.showlegend = true;
    h.layout.legend = struct('x', 0.1, 'y', 0.95);

    h.PlotOptions.SaveFolder = outdir;
    h.PlotOptions.FileName = 'eye_events_traffic_decision';
    plotlyoffline(h);

    if show_plots
        web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
    end

    % Correct/incorrect
    h = plot_timelocked(params, summary.stats.traffic_decision.correct, {'Incorrect','Correct'}, ...
                        plotly_layouts.layouts.boxplot, sig_dims, sig_clrs);
    h.layout.title = 'PD Locked to Traffic Decision: Correct v. Incorrect';
    h.layout.yaxis.range = params.eye.events.traffic_decision.plots.ylims; % [-0.5 2]; % [min(y(:)) max(y(:))]*1.3;
    h.layout.showlegend = true;
    h.layout.legend = struct('x', 0.1, 'y', 0.95);

    h.PlotOptions.SaveFolder = outdir;
    h.PlotOptions.FileName = 'eye_events_traffic_decision_correct';
    plotlyoffline(h);

    if show_plots
        web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
    end
    
    % High/Low confidence
    h = plot_timelocked(params, summary.stats.traffic_decision.confidence, {'Low','High'}, ...
                        plotly_layouts.layouts.boxplot, sig_dims, sig_clrs);
    h.layout.title = 'PD Locked to Traffic Decision: Confidence';
    h.layout.yaxis.range = params.eye.events.traffic_decision.plots.ylims; % [-0.5 2]; % [min(y(:)) max(y(:))]*1.3;
    h.layout.showlegend = true;
    h.layout.legend = struct('x', 0.1, 'y', 0.95);

    h.PlotOptions.SaveFolder = outdir;
    h.PlotOptions.FileName = 'eye_events_traffic_decision_confidence';
    plotlyoffline(h);

    if show_plots
        web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
    end
    
    % 1-/2-back
    h = plot_timelocked(params, summary.stats.traffic_decision.order, {'Two-Back','One-Back'}, ...
                        plotly_layouts.layouts.boxplot, sig_dims, sig_clrs);
    h.layout.title = 'PD Locked to Traffic Decision: One- vs. Two-Back';
    h.layout.yaxis.range = params.eye.events.traffic_decision.plots.ylims; % [-0.5 2]; % [min(y(:)) max(y(:))]*1.3;
    h.layout.showlegend = true;
    h.layout.legend = struct('x', 0.1, 'y', 0.95);

    h.PlotOptions.SaveFolder = outdir;
    h.PlotOptions.FileName = 'eye_events_traffic_decision_order';
    plotlyoffline(h);

    if show_plots
        web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
    end
    
    % High/Low confidence
    h = plot_timelocked(params, summary.stats.traffic_decision.confidence, {'Low','High'}, ...
                        plotly_layouts.layouts.boxplot, sig_dims, sig_clrs);
    h.layout.title = 'PD Locked to Traffic Decision: Confidence';
    h.layout.yaxis.range = params.eye.events.traffic_decision.plots.ylims; % [-0.5 2]; % [min(y(:)) max(y(:))]*1.3;
    h.layout.showlegend = true;
    h.layout.legend = struct('x', 0.1, 'y', 0.95);

    h.PlotOptions.SaveFolder = outdir;
    h.PlotOptions.FileName = 'eye_events_traffic_decision_confidence';
    plotlyoffline(h);

    if show_plots
        web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
    end
    
end

% Regression analyses
if ~isempty(params.eye.events.covariates.glms)
    
     for i = 1 : length(params.eye.events.covariates.glms)
        
        glm_i = params.eye.events.covariates.glms{i};
        event_types = fieldnames(summary.stats.covariate_glms.(glm_i.name));
        
        for j = 1 : length(event_types)
            event_type = event_types{j};
            summary_j = summary.stats.covariate_glms.(glm_i.name).(event_type);
            
            % Time series
            h = plot_timelocked_stats(params, summary_j, ...
                                      {''}, ...
                                      plotly_layouts.layouts.boxplot, sig_dims, sig_clrs);

            h.layout.title = sprintf('GLM: %s - %s', glm_i.name, summary_j.event_title);
            h.layout.yaxis.range = params.eye.events.covariates.plots.ylims; % [-0.5 2]; % [min(y(:)) max(y(:))]*1.3;
%             h.layout.showlegend = true;
%             h.layout.legend = struct('x', 0.1, 'y', 3.0);
            
            h.PlotOptions.SaveFolder = outdir;
            h.PlotOptions.FileName = sprintf('eye_events_glm_%s_%s', glm_i.name, event_type);
            plotlyoffline(h);

            if show_plots
                web(sprintf('file://%s/%s.html', outdir, h.PlotOptions.FileName), '-new', '-notoolbar');
            end




            % Scatterplots
    
        
        end
    
     end
    
    
end




end