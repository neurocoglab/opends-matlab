alpha_level = 0.7;
idx_channel = 36;

mask = squeeze(double(stats.left_change.timelockstats.mask(idx_channel,:,:)));
mask(mask==0) = 0.7;
fmask = flipud(mask);
stat = squeeze(stats.left_change.timelockstats.stat(idx_channel,:,:));
time = stats.left_change.timelockstats.time;
freq = stats.left_change.timelockstats.freq;

figure;
pcolor(time, freq, stat);
shading flat;


figure;
pcolor(time, freq, stat);

hold on;
h=pcolor(time, freq, mask);
h.AlphaData=mask*200;
shading flat;


figure;
h=pcolor(time, freq, stat);
h.AlphaData=mask;
shading flat;

figure;
h = imagesc(flipud(stat));
h.AlphaData = fmask;
hold on;
[~,h]=imcontour(fmask,[0,alpha_level,1]);
h.Color = 'k';

figure;
h = imagesc(flipud(stat));
h.AlphaData = fmask;
hold on;

bmask = flipud(squeeze(stats.left_change.timelockstats.mask(idx_channel,:,:)));
[i,j] = find(bmask);
b = contour(i,j);
figure,plot(i(b),j(b));

figure;
h = imagesc(flipud(stat));
h.AlphaData = fmask;
hold on;

bb = bwboundaries(bmask);
for i = 1 : length(bb)
   plot(bb{i}(:,2),bb{i}(:,1),'k');
   hold on;
end

figure;
bmask = flipud(squeeze(stats.left_change.timelockstats.mask(idx_channel,:,:)));
h = imagesc(flipud(stat));
h.AlphaData = fmask;
hold on;
bw = bwperimtrace(bmask');
for i = 1 : length(bw)
   plot(bw{i}(:,1),bw{i}(:,2),'Color',[130 38 31]/255);
   hold on;
end

figure,imagesc(bmask);

plot_params = [];
plot_params.xtick = -8:1:6;
plot_params.ytick = [[2:2:22],[24:4:40]];
plot_params.title = 'Overtake onset';
plot_params.channels = {'Fz'};
plot_params.transparency = 0.7;
plot_params.window_size = [800 600];
plot_params.clim = [-6 6];
plot_params.colorbar = true;

h = plot_timefreq_img(stats.left_change.timelockstats, plot_params);

plot_params.xtick = -5:1:5;
plot_params.ytick = [[2:2:22],[24:4:40]];
plot_params.title = 'Overtake offset';
h = plot_timefreq_img(stats.right_change.timelockstats, plot_params);


