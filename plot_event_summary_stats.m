function [] = plot_event_summary_stats ( summary, params )

% Overtake
ci_colour = [0.1 0.1 0.9];
zero_colour = [1 0 0];

h = figure;
set(h,'Color','w');

overtake = summary.overtake.pupil;
pdm = mean(overtake,1)';
pdv = std(overtake,0,1)';

plot_ci_filled(results.overtake.t', pdm, pdm-pdv, pdm+pdv, ci_colour); 

hold on;

hh = line([0 0],[-1000 1000]);
set(hh,'Color',zero_colour);
set(hh,'LineWidth',1.5);

ylim([-2 2]);

xlabel('Time from overtake event (s)');
ylabel('Relative pupil diameter (z-score)');

% Lane change left
ci_colour = [0.1 0.1 0.9];
zero_colour = [1 0 0];

h = figure;
set(h,'Color','w');

overtake = summary.left_change.pupil;
pdm = mean(overtake,1)';
pdv = std(overtake,0,1)';

plot_ci_filled(results.left_change.t', pdm, pdm-pdv, pdm+pdv, ci_colour); 

hold on;

hh = line([0 0],[-1000 1000]);
set(hh,'Color',zero_colour);
set(hh,'LineWidth',1.5);

ylim([-2 2]);

xlabel('Time from left lane change (s)');
ylabel('Relative pupil diameter (z-score)');

% Lane change right
ci_colour = [0.1 0.1 0.9];
zero_colour = [1 0 0];

h = figure;
set(h,'Color','w');

overtake = summary.right_change.pupil;
pdm = mean(overtake,1)';
pdv = std(overtake,0,1)';

plot_ci_filled(results.right_change.t', pdm, pdm-pdv, pdm+pdv, ci_colour); 

hold on;

hh = line([0 0],[-1000 1000]);
set(hh,'Color',zero_colour);
set(hh,'LineWidth',1.5);

ylim([-2 2]);

xlabel('Time from right lane change (s)');
ylabel('Relative pupil diameter (z-score)');


end