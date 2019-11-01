%% Simple script to plot correlated values

load('drive_game_n40.mat');

h = figure;
set(h,'Color','w');

hs = scatter(sim_score,drive_game_exp,100,'b','filled');

set(gca,'FontSize',11);

hl = lsline;
set(hl,'LineWidth',3);
set(hl,'Color','b');
corr(sim_score,drive_game_exp);
hx = xlabel('Simulation Score');
hy = ylabel('Driving/Gaming Experience');

set(hx,'FontSize',20);
set(hy,'FontSize',20);

%% Plots time v. luminance

h = figure;

set(h,'Color','w');

hh = plot(Time,Luminance,'b');
ylim([0 1]);
xlim([0 10]);

set(hh,'LineWidth',2);

