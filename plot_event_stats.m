function [] = plot_event_stats ( results, params )

% Overtake
ci_colour = [0.1 0.1 0.9];
zero_colour = [1 0 0];

h = figure;
set(h,'Color','w');

overtake = results.overtake.tlocked_bl;
pdm = mean(overtake,1,'omitnan')';
pdv = std(overtake,0,1,'omitnan')';
pderr = pdv / results.N;

plot_ci_filled(results.overtake.t', pdm, pdm-pderr, pdm+pderr, ci_colour);
set(gca, 'FontName', 'Arial Narrow');
set(gca, 'FontSize', params.events.plots.tick_label_size)

hold on;

hh = line([0 0],[-1000 1000]);
set(hh,'Color',zero_colour);
set(hh,'LineWidth',1.5);

if ~isempty(params.events.plots.xlim)
    xlim(params.events.plots.xlim);
end

if ~isempty(params.events.plots.ylim)
    ylim(params.events.plots.ylim);
else
    buffer = range(pdm)*0.2;
    ylim([min(pdm)-buffer max(pdm)+buffer]);
end

% ylim(params.events.overtake.ylim); %[-2 2]);
if params.events.plots.label_axes
    h=xlabel('Time from overtake event (s)');
    set(h,'FontSize',params.events.plots.label_size);
    h=ylabel('Relative pupil diameter (z-score)');
    set(h,'FontSize',params.events.plots.label_size);
end

if params.events.save_plots
    print(sprintf('%s/%sovertake_tlocked.png',params.output_dir, ...
                                              params.events.plot_prefix), ...
          '-dpng');
end


% Overtake by difficulty
ci_clrs=[[0.1 0.1 0.9]; ...
         [0.5 0.5 0.1]];

h = figure;
set(h,'Color','w');

for i = 1 : 2
    overtake = results.overtake.tlocked_bl(results.overtake.diffs==i,:);
    pdm = mean(overtake,1,'omitnan')';
    pdv = std(overtake,0,1,'omitnan')';
    pderr = pdv / results.N;
    
    plot_ci_filled(results.overtake.t', pdm, pdm-pderr, pdm+pderr, ci_clrs(i,:));
    hold on;
end

hold on;

hh = line([0 0],[-1000 1000]);
set(hh,'Color',zero_colour);
set(hh,'LineWidth',1.5);

if ~isempty(params.events.plots.xlim)
    xlim(params.events.plots.xlim);
end

if ~isempty(params.events.plots.ylim)
    ylim(params.events.plots.ylim);
else
    buffer = range(pdm)*0.2;
    ylim([min(pdm)-buffer max(pdm)+buffer]);
end

% ylim(params.events.overtake.ylim); %[-2 2]);
if params.events.plots.label_axes
    h=xlabel('Time from overtake event (s)');
    set(h,'FontSize',params.events.plots.label_size);
    h=ylabel('Relative pupil diameter (z-score)');
    set(h,'FontSize',params.events.plots.label_size);
end

% Show legend
l1 = plot([NaN,NaN], 'color', ci_clrs(1,:));
l2 = plot([NaN,NaN], 'color', ci_clrs(2,:));
legend([l1, l2], {'Easy', 'Difficult'});

if params.events.save_plots
    print(sprintf('%s/%sovertake_tlocked_diff.png',params.output_dir, ...
                                              params.events.plot_prefix), ...
          '-dpng');
end


% Lane change left
ci_colour = [0.1 0.1 0.9];
zero_colour = [1 0 0];

h = figure;
set(h,'Color','w');

overtake = results.left_change.tlocked_bl;
pdm = mean(overtake,1,'omitnan')';
pdv = std(overtake,0,1,'omitnan')';
pderr = pdv / results.N;

plot_ci_filled(results.left_change.t', pdm, pdm-pderr, pdm+pderr, ci_colour); 
set(gca, 'FontName', 'Arial Narrow');
set(gca, 'FontSize', params.events.plots.tick_label_size)

hold on;

hh = line([0 0],[-1000 1000]);
set(hh,'Color',zero_colour);
set(hh,'LineWidth',1.5);

if ~isempty(params.events.plots.xlim)
    xlim(params.events.plots.xlim);
end

if ~isempty(params.events.plots.ylim)
    ylim(params.events.plots.ylim);
else
    buffer = range(pdm)*0.2;
    ylim([min(pdm)-buffer max(pdm)+buffer]);
end
% ylim(params.events.left_change.ylim); %[-2 2]);

if params.events.plots.label_axes
    h=xlabel('Time from left lane change (s)');
    set(h,'FontSize',params.events.plots.label_size);
    h=ylabel('Relative pupil diameter (z-score)');
    set(h,'FontSize',params.events.plots.label_size);
end

if params.events.save_plots
    print(sprintf('%s/%sleft_change_tlocked.png',params.output_dir, ...
                                                 params.events.plot_prefix), ...
          '-dpng');
end


% Lane change left by difficulty
ci_clrs=[[0.1 0.1 0.9]; ...
         [0.5 0.5 0.1]];

h = figure;
set(h,'Color','w');

for i = 1 : 2
    overtake = results.left_change.tlocked_bl(results.left_change.diffs==i,:);
    pdm = mean(overtake,1,'omitnan')';
    pdv = std(overtake,0,1,'omitnan')';
    pderr = pdv / results.N;
    
    plot_ci_filled(results.left_change.t', pdm, pdm-pderr, pdm+pderr, ci_clrs(i,:));
    hold on;
end

hold on;

hh = line([0 0],[-1000 1000]);
set(hh,'Color',zero_colour);
set(hh,'LineWidth',1.5);

if ~isempty(params.events.plots.xlim)
    xlim(params.events.plots.xlim);
end

if ~isempty(params.events.plots.ylim)
    ylim(params.events.plots.ylim);
else
    buffer = range(pdm)*0.2;
    ylim([min(pdm)-buffer max(pdm)+buffer]);
end

% ylim(params.events.overtake.ylim); %[-2 2]);
if params.events.plots.label_axes
    h=xlabel('Time from left lane change (s)');
    set(h,'FontSize',params.events.plots.label_size);
    h=ylabel('Relative pupil diameter (z-score)');
    set(h,'FontSize',params.events.plots.label_size);
end

% Show legend
l1 = plot([NaN,NaN], 'color', ci_clrs(1,:));
l2 = plot([NaN,NaN], 'color', ci_clrs(2,:));
legend([l1, l2], {'Easy', 'Difficult'});

if params.events.save_plots
    print(sprintf('%s/%sleft_change_tlocked_diff.png',params.output_dir, ...
                                              params.events.plot_prefix), ...
          '-dpng');
end



% Lane change right
ci_colour = [0.1 0.1 0.9];
zero_colour = [1 0 0];

h = figure;
set(h,'Color','w');

overtake = results.right_change.tlocked_bl;
pdm = mean(overtake,1,'omitnan')';
pdv = std(overtake,0,1,'omitnan')';
pderr = pdv / results.N;

plot_ci_filled(results.right_change.t', pdm, pdm-pderr, pdm+pderr, ci_colour); 
set(gca, 'FontName', 'Arial Narrow');
set(gca, 'FontSize', params.events.plots.tick_label_size)

hold on;

hh = line([0 0],[-1000 1000]);
set(hh,'Color',zero_colour);
set(hh,'LineWidth',1.5);

if ~isempty(params.events.plots.xlim)
    xlim(params.events.plots.xlim);
end

if ~isempty(params.events.plots.ylim)
    ylim(params.events.plots.ylim);
else
    buffer = range(pdm)*0.2;
    ylim([min(pdm)-buffer max(pdm)+buffer]);
end
% ylim(params.events.right_change.ylim); %[-2 2]);

if params.events.plots.label_axes
    h=xlabel('Time from right lane change (s)');
    set(h,'FontSize',params.events.plots.label_size);
    h=ylabel('Relative pupil diameter (z-score)');
    set(h,'FontSize',params.events.plots.label_size);
end

if params.events.save_plots
    print(sprintf('%s/%sright_change_tlocked.png',params.output_dir, ...
                                                  params.events.plot_prefix), ...
          '-dpng');
end


% Lane change left
ci_colour = [0.1 0.1 0.9];
zero_colour = [1 0 0];

h = figure;
set(h,'Color','w');

overtake = results.left_change.tlocked_bl;
pdm = mean(overtake,1,'omitnan')';
pdv = std(overtake,0,1,'omitnan')';
pderr = pdv / results.N;

plot_ci_filled(results.left_change.t', pdm, pdm-pderr, pdm+pderr, ci_colour); 
set(gca, 'FontName', 'Arial Narrow');
set(gca, 'FontSize', params.events.plots.tick_label_size)

hold on;

hh = line([0 0],[-1000 1000]);
set(hh,'Color',zero_colour);
set(hh,'LineWidth',1.5);

if ~isempty(params.events.plots.xlim)
    xlim(params.events.plots.xlim);
end

if ~isempty(params.events.plots.ylim)
    ylim(params.events.plots.ylim);
else
    buffer = range(pdm)*0.2;
    ylim([min(pdm)-buffer max(pdm)+buffer]);
end
% ylim(params.events.left_change.ylim); %[-2 2]);

if params.events.plots.label_axes
    h=xlabel('Time from left lane change (s)');
    set(h,'FontSize',params.events.plots.label_size);
    h=ylabel('Relative pupil diameter (z-score)');
    set(h,'FontSize',params.events.plots.label_size);
end

if params.events.save_plots
    print(sprintf('%s/%sleft_change_tlocked.png',params.output_dir, ...
                                                 params.events.plot_prefix), ...
          '-dpng');
end


% Lane change right by difficulty
ci_clrs=[[0.1 0.1 0.9]; ...
         [0.5 0.5 0.1]];

h = figure;
set(h,'Color','w');

for i = 1 : 2
    overtake = results.right_change.tlocked_bl(results.right_change.diffs==i,:);
    pdm = mean(overtake,1,'omitnan')';
    pdv = std(overtake,0,1,'omitnan')';
    pderr = pdv / results.N;
    
    plot_ci_filled(results.right_change.t', pdm, pdm-pderr, pdm+pderr, ci_clrs(i,:));
    hold on;
end

hold on;

hh = line([0 0],[-1000 1000]);
set(hh,'Color',zero_colour);
set(hh,'LineWidth',1.5);

if ~isempty(params.events.plots.xlim)
    xlim(params.events.plots.xlim);
end

if ~isempty(params.events.plots.ylim)
    ylim(params.events.plots.ylim);
else
    buffer = range(pdm)*0.2;
    ylim([min(pdm)-buffer max(pdm)+buffer]);
end

% ylim(params.events.overtake.ylim); %[-2 2]);
if params.events.plots.label_axes
    h=xlabel('Time from right lane change (s)');
    set(h,'FontSize',params.events.plots.label_size);
    h=ylabel('Relative pupil diameter (z-score)');
    set(h,'FontSize',params.events.plots.label_size);
end

% Show legend
l1 = plot([NaN,NaN], 'color', ci_clrs(1,:));
l2 = plot([NaN,NaN], 'color', ci_clrs(2,:));
legend([l1, l2], {'Easy', 'Difficult'});

if params.events.save_plots
    print(sprintf('%s/%sright_change_tlocked_diff.png',params.output_dir, ...
                                              params.events.plot_prefix), ...
          '-dpng');
end



% Saccades
ci_colour = [0.1 0.1 0.9];
zero_colour = [1 0 0];

h = figure;
set(h,'Color','w');

overtake = results.saccades.tlocked;
pdm = mean(overtake,1,'omitnan')';
pdv = std(overtake,0,1,'omitnan')';
pderr = pdv / results.N;

plot_ci_filled(results.saccades.t', pdm, pdm-pderr, pdm+pderr, ci_colour); 
set(gca, 'FontName', 'Arial Narrow');
set(gca, 'FontSize', params.events.plots.tick_label_size)

hold on;

hh = line([0 0],[-1000 1000]);
set(hh,'Color',zero_colour);
set(hh,'LineWidth',1.5);

if ~isempty(params.events.plots.xlim)
    xlim(params.events.plots.xlim);
end

if ~isempty(params.events.plots.ylim)
    ylim(params.events.plots.ylim);
else
    buffer = range(pdm)*0.2;
    ylim([min(pdm)-buffer max(pdm)+buffer]);
end
% ylim(params.events.saccades.ylim); %);

if params.events.plots.label_axes
    h=xlabel('Saccade middle (s)');
    set(h,'FontSize',params.events.plots.label_size);
    h=ylabel('Relative pupil diameter (z-score)');
    set(h,'FontSize',params.events.plots.label_size);
end

if params.events.save_plots
    print(sprintf('%s/%ssaccades_tlocked.png',params.output_dir, ...
                                              params.events.plot_prefix), ...
          '-dpng');
end

% Blinks
ci_colour = [0.1 0.1 0.9];
zero_colour = [1 0 0];

h = figure;
set(h,'Color','w');

overtake = results.blinks.tlocked;
pdm = mean(overtake,1,'omitnan')';
pdv = std(overtake,0,1,'omitnan')';
pderr = pdv / results.N;

plot_ci_filled(results.blinks.t', pdm, pdm-pderr, pdm+pderr, ci_colour); 
set(gca, 'FontName', 'Arial Narrow');
set(gca, 'FontSize', params.events.plots.tick_label_size)

hold on;

hh = line([0 0],[-1000 1000]);
set(hh,'Color',zero_colour);
set(hh,'LineWidth',1.5);

if ~isempty(params.events.plots.xlim)
    xlim(params.events.plots.xlim);
end

if ~isempty(params.events.plots.ylim)
    ylim(params.events.plots.ylim);
else
    buffer = range(pdm)*0.2;
    ylim([min(pdm)-buffer max(pdm)+buffer]);
end
% ylim(params.events.blinks.ylim); %);

if params.events.plots.label_axes
    h=xlabel('Blink offset (s)');
    set(h,'FontSize',params.events.plots.label_size);
    h=ylabel('Relative pupil diameter (z-score)');
    set(h,'FontSize',params.events.plots.label_size);
end

if params.events.save_plots
    print(sprintf('%s/%sblinks_tlocked.png',params.output_dir, ...
                                            params.events.plot_prefix), ...
          '-dpng');
end

% Random
if isfield(results, 'perms')
    ci_colour = [0.1 0.1 0.9];
    zero_colour = [1 0 0];

    h = figure;
    set(h,'Color','w');

    overtake = results.perms.tlocked;
    pdm = mean(overtake,1,'omitnan')';
    pdv = std(overtake,0,1,'omitnan')';
    pderr = pdv / results.N;

    plot_ci_filled(results.perms.t', pdm, pdm-pderr, pdm+pderr, ci_colour); 
    set(gca, 'FontName', 'Arial Narrow');
    set(gca, 'FontSize', params.events.plots.tick_label_size)

    hold on;

    hh = line([0 0],[-1000 1000]);
    set(hh,'Color',zero_colour);
    set(hh,'LineWidth',1.5);

    if ~isempty(params.events.plots.xlim)
        xlim(params.events.plots.xlim);
    end

    if ~isempty(params.events.plots.ylim)
        ylim(params.events.plots.ylim);
    else
        buffer = range(pdm)*0.2;
        ylim([min(pdm)-buffer max(pdm)+buffer]);
    end
    % ylim(params.events.perms.ylim); %);

    if params.events.plots.label_axes
        h = xlabel('Random index offset (s)');
        set(h,'FontSize',params.events.plots.label_size);
        h = ylabel('Relative pupil diameter (z-score)');
        set(h,'FontSize',params.events.plots.label_size);
    end

    if params.events.save_plots
        print(sprintf('%s/%srandom_tlocked.png',params.output_dir, ...
                                                  params.events.plot_prefix), ...
              '-dpng');
    end
end

end