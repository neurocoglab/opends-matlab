function [ h ] = plot_artifacts_eeg(params, data, out2file)

% Plot EEG and blinks together
% Time variables must have already been aligned

if nargin < 3
    out2file = false;
end

channels = params.eeg.artifacts.plots.channels;
idx_channels = [];
for i = 1 : length(channels)
    idx_channels = [idx_channels find(strcmp(data.eeg.ft.label,channels{i}))];
end

if out2file
    h = figure('visible','off');
else
    h = figure;
end

outdir = sprintf( '%s/%s', params.io.output_dir, data.subject );
results_file = sprintf('%s/results_preproc_eye.mat', outdir);
T = load( results_file, 'data' );
data.sim = T.data.sim;
clear T;

h.Color = [1 1 1];
x_eye = seconds(data.eye.t / 1000);

N_eye = length(data.eye.t);
idx_keep = 1 : N_eye;

t_end = data.sim.simended.values.Millis;
if ~isnan(t_end)
    t_end = seconds((t_end-data.sim.t_start) / 1000);
    idx_keep = find(x_eye <= t_end);
    x_eye = x_eye(idx_keep);
end

idx_sacc = data.eye.saccades.saccades(:,1:2);
idx_sacc = idx_sacc(~any(idx_sacc>max(idx_keep),2),:);
for i = 1 : size(idx_sacc,1)
    xs = x_eye(idx_sacc(i,:));
    hh = fill([xs(1) xs(1) xs(2) xs(2)],[-1000 1000 1000 -1000],[0.5 0.1 0.1]);     
    hh.EdgeColor = 'none';
    alpha(hh, 0.3);
    hold on;
end

intervals = data.eye.blinks.intervals;
intervals = intervals(~any(intervals>max(idx_keep),2),:);

for i = 1 : length(intervals)
    xs = [x_eye(intervals(i,1)) x_eye(intervals(i,2))];
    hh = fill([xs(1) xs(1) xs(2) xs(2)],[-1000 1000 1000 -1000],[0.1 0.1 0.5]);     
    hh.EdgeColor = 'none';
    alpha(hh, 0.3);
    hold on;
end

intervals = data.eye.blinks.intervals;
intervals = intervals(~any(intervals>max(idx_keep),2),:);

for i = 1 : length(intervals)
    xs = [x_eye(intervals(i,1)) x_eye(intervals(i,2))];
    hh = fill([xs(1) xs(1) xs(2) xs(2)],[-1000 1000 1000 -1000],[0.1 0.1 0.5]);     
    hh.EdgeColor = 'none';
    alpha(hh, 0.3);
    hold on;
end


L = 0;
lens = zeros(length(data.eeg.ft.time),1);
for c = 1 : length(data.eeg.ft.time)
    lens(c) = length(data.eeg.ft.time{c});
    L = L + lens(c) + 1;
end
x_eeg = duration(nan(L,3));
y_eeg = zeros(length(idx_channels),L);
idx = 1;
for c = 1 : length(data.eeg.ft.time)
    x_eeg(idx:idx+lens(c)-1) = seconds(data.eeg.ft.time{c});
    x_eeg(idx+lens(c)) = x_eeg(idx+lens(c)-1) + seconds(0.00001);
    yy = data.eeg.ft.trial{c};
    if ~isempty(idx_channels)
       y_eeg(:,idx:idx+lens(c)-1) = yy(idx_channels,:); 
    end
    y_eeg(:,idx+lens(c)) = nan;
    idx = idx+lens(c)+1;
end

idx_keep2 = x_eeg >= x_eye(1) & x_eeg <= x_eye(end);
x_eeg = x_eeg(idx_keep2);
y_eeg = y_eeg(:,idx_keep2);

N_ch = size(y_eeg,1);
itr = 300 / N_ch;
gain = itr / params.eeg.artifacts.plots.stdev; % Std devs per lane

for i = 1 : N_ch
   y_eeg(i,:) = zscore_nan(y_eeg(i,:)) * gain - (i-1)*itr;
end

hh = plot(x_eeg, y_eeg);
ypos = zeros(N_ch,1);
ystr = cell(N_ch,1);
for i = 1 : N_ch
   ypos(i) = nanmean(y_eeg(i,:));
   ystr(i) = data.eeg.ft.label(idx_channels(i));
end
ax = gca;

pd = data.eye.blinks.diam;
y_posx = 150 + zscore(pd(idx_keep)) * 10;
hh = plot(x_eye, y_posx, 'b', 'LineWidth', 2);

y_posx = 110 + zscore(data.eye.pos_x(idx_keep)) * 10;
hh = plot(x_eye, y_posx, 'Color', [0 .3 0], 'LineWidth', 2);

y_posy = 75 + zscore(data.eye.pos_y(idx_keep)) * 10;
hh = plot(x_eye, y_posy, 'm', 'LineWidth', 2);

xlim(seconds([0 10]));
ylim([-300 200]);

xlabel('Time (mm:ss)');
ax.XAxis.Label.FontSize = 20;

ax.YTick=flip(ypos);
ax.YTickLabel=flip(ystr);
ax.XAxis.TickLabelFormat = 'mm:ss';

ax.FontSize=17;

if out2file
    fig_file = sprintf('%s/figures/artifacts_eeg', outdir);
    saveas(h,sprintf('%s.fig', fig_file));
    saveas(h,sprintf('%s.png', fig_file));
    close(h);
end

    
function Z = zscore_nan(X)
        
    Z = (X-nanmean(X)) ./ nanstd(X);

end


end