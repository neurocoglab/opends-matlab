function [ h ] = plot_epoch_stats ( results )

% Full session
h = figure;

M = [results.baseline.pupil;results.nobaseline.pupil;results.passing.pupil];
grp = [zeros(length(results.baseline.pupil),1);ones(length(results.nobaseline.pupil),1); ...
            2*ones(length(results.passing.pupil),1)];

boxplot(M,grp,'outliersize',3,'symbol','');

% ylim([-5 5]);

hh = title('Pupil diameter - Full session');
set (hh,'FontSize',16);
       
set(gca, 'XTickLabel', [{'Baseline'},{'Non-baseline'},{'Passing'}]);
set (gca,'FontSize',14);

hh = ylabel('Pupil diameter (z-score)');
set (hh,'FontSize',14);

% By difficulty
h = figure;
M = results.passing.pupil;
grp =results.passing.diff;
boxplot(M,grp,'outliersize',3,'symbol','');
hh = title('Pupil diameter - Passing difficulty');
set (hh,'FontSize',16);
set(gca, 'XTickLabel', [{'Easy'},{'Moderate'},{'Difficult'}]);
set (gca,'FontSize',14);
hh = ylabel('Pupil diameter (z-score)');
set (hh,'FontSize',14);


% Baselines over time
h = figure;

N = length(results.intervals.baseline.pupil);
x = results.intervals.baseline.times / 60000;
y = zeros(N,1);
err = zeros(N,1);

for i = 1 : N
    pdi = results.intervals.baseline.pupil{i};
    y(i) = mean(pdi);
    err(i) = std(pdi);  
end
    
colours = [[0 0 1];[1 0 0]];
lighten = 0.8;

hp = zeros(2, 1);
hp(1) = plot_ci_filled( x, y, y+err, y-err, colours(1,:), lighten );
%hh = errorbar(x,y,err,'b');
hold on;

N = length(results.intervals.passing.pupil);
x = results.intervals.passing.times / 60000;
y = zeros(N,1);
err = zeros(N,1);

for i = 1 : N
    pdi = results.intervals.passing.pupil{i};
    y(i) = mean(pdi);
    err(i) = std(pdi);  
end

%errorbar(x,y,err,'r');
hp(2) = plot_ci_filled( x, y, y+err, y-err, colours(2,:), lighten );

hh = legend(hp,'Baseline','Passing');
set (hh,'FontSize',14);

% 
hh = title('Pupil diameter - Baseline Intervals');
set (hh,'FontSize',16);
%        
hh = ylabel('Pupil diameter (z-score)');
set (hh,'FontSize',14);

hh = xlabel('Time (min)');
set (hh,'FontSize',14);

% By cycle/round






end