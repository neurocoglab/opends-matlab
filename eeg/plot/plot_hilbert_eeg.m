function [ h ] = plot_hilbert_eeg ( params, data, results, out2file )

% PLOT_HILBERT_EEG Plots the results of Hilbert processing as a set of
%   time series, plotted in rows, showing pupil diameter and each Hilbert
%   envelope over time
%
% params:       Parameter struct
% results:      Results struct (must include hilbert and eye/sim results
% out2file:     If true, saves the plot to a file rather than showing it in
%               a new window
%

% Interval colours
gap_clr = [params.eye.gaps.plots.color 0.5];
blink_clr = [params.eye.blinks.plots.color 0.5];
pass_clr = [params.sim.plots.pass_color 0.5];

% Which bands to plot?
bands = params.eeg.hilbert.plot.bands;
if strcmp(bands, 'all')
   bands = results.eeg.hilbert.bands; 
end
N_bands = height(bands);

% Which channels to average?
channels = params.eeg.hilbert.plot.channels;
if strcmp(channels, 'all')
   channels = results.eeg.hilbert.channels; 
end
N_channels = length(channels);

if N_channels == length(results.eeg.hilbert.channels)
   idx_channels = 1 :  N_channels;
else
   idx_channels = zeros(N_channels,1);
   for i = 1 : N_channels
       idx_channels(i) = find(strcmp(results.eeg.hilbert.channels, channels{i}),1);
   end
end

% Line colours
color_pd = [0 0 1];
colors = params.general.plots.plotly_colors;

% Line plot parameters
scale_pd = params.eeg.hilbert.plot.scale_pupil;
scale_env = params.eeg.hilbert.plot.scale_envelope;
scale_g = 1 / params.eeg.hilbert.plot.scale_general;

pdz = data.eye.diam;
if isfield(data.eye, 'blinks')
   pdz = data.eye.blinks.diam; 
end

pdz = zscore(pdz);

% Smoothing?
if params.eeg.hilbert.plot.smooth_pupil > 0
   pdz = smooth(pdz, params.eeg.hilbert.plot.smooth_pupil); 
end


t_eye = data.eye.t / 60000;
% t_eeg = results.eeg.hilbert.time / 60;

Ylim_pd = 2 * scale_pd;
Ylim_env = 2 * scale_env * N_bands;
Ylim_total = scale_g / 2 * [-Ylim_pd - Ylim_env, Ylim_pd];

Ysep_pd = scale_g * Ylim_pd / 2;
Ysep_env = scale_g * Ylim_env / N_bands / 2;

if out2file
    h = figure('visible','off');
else
    h = figure;
end

% Plot data gaps
for j = 1 : size(data.eye.tgap,1)
   x1 = t_eye(data.eye.tgap(j,1));
   x2 = x1 + data.eye.tgap(j,2) / 60000;
   if x2 > x1
        a = rectangle('Position',[x1 -500 x2-x1 1000], ...
             'EdgeColor','w', ...
             'FaceColor', gap_clr); 
   end
end

hold on;

% Blinks
for j = 1 : size(data.eye.blinks.intervals,1)
   x1 = t_eye(data.eye.blinks.intervals(j,1));
   try
   x2 = t_eye(data.eye.blinks.intervals(j,2));
   catch
       a=0;
   end
   a = rectangle('Position',[x1 -500 x2-x1 1000], ...
                  'EdgeColor','w', ...
                  'FaceColor', blink_clr); 
end

% Plot lane changes
max_interval = 30000; % Half a minute

left = [data.sim.sim2track.left_change_times, ...
        true(length(data.sim.sim2track.left_change_times),1)];
right = [data.sim.sim2track.right_change_times, ...
         false(length(data.sim.sim2track.right_change_times),1)];

left_right = [left;right];
[~,idx] = sort(left_right(:,1));
left_right = left_right(idx,:);

this_left = -1;
this_right = -1;
for j = 1 : length(left_right)
   if left_right(j,2)
      % Is change to left 
      this_left = left_right(j,1);
   else
      % Is change to right
      this_right = left_right(j,1);
      if this_left > 0 && this_right-this_left < max_interval
          % Valid passing segment, plot
          x1 = this_left / 60000;
          x2 = this_right / 60000;
          if x2 > x1
          a = rectangle('Position',[x1 -500 x2-x1 1000], ...
                     'EdgeColor','w', ...
                     'FaceColor', pass_clr); 
          end
      end
      % Reset
      this_left = -1;
      this_right = -1;
   end
end

yticks = [0];
ylabels = {};

% Plot envelopes first
for i = N_bands : -1 : 1
   
    band = bands.Band{i};
    idx_i = find(strcmp(results.eeg.hilbert.bands.Band, band));
    
    envelopes_i = results.eeg.hilbert.envelopes{idx_i};
    y = zeros(length(idx_channels),length(envelopes_i{1}));
    for j = 1 : length(idx_channels)
       y(j,:) = abs(envelopes_i{idx_channels(j)});
    end

    t_eeg = results.eeg.hilbert.time{1}{1};
    y = mean(y,1); y = y(:);
    y(y>0) = zscore(y(y>0));
    y(y==0) = nan;
    y_offset = Ysep_pd + (i-1) * Ysep_env;
    y = y - y_offset;
    yticks = [yticks -y_offset];
    
    hh = plot(t_eeg, y*scale_env);
    hold on;
    clr = colors{i};
    if clr(1)=='#'
       clr = sscanf(clr(2:end),'%2x%2x%2x',[1 3])/255;
    end
    hh.Color = clr;
    
    ylabel = sprintf('\\color[rgb]{%1.2f,%1.2f,%1.2f} %s', clr(1), clr(2), clr(3), band);
    ylabels = [ylabels {ylabel}];

end

% Plot PD
hh = plot(t_eye, pdz*scale_pd);
clr = color_pd;
if clr(1)=='#'
   clr = sscanf(clr(2:end),'%2x%2x%2x',[1 3])/255;
end
ylabel = sprintf('\\color[rgb]{%1.2f,%1.2f,%1.2f} PD', clr(1), clr(2), clr(3));
ylabels = [{ylabel} ylabels];
hh.Color = clr;
hh.LineWidth = 2;

xlabel('Time (mins)');
xlim(params.eeg.hilbert.plot.xlims);

h.Position = [200 200 1500 750];
h.Color = 'w';

hh = title(sprintf('Hilbert envelopes for %s', results.subject));
hh.FontSize = 18;
hh.FontWeight = 'bold';

ylim(Ylim_total);

ax = gca;

[sy,idx] = sort(yticks);
ax.YTick = sy;
ax.YTickLabel = ylabels(idx);
ax.YAxis.FontSize = 16;
ax.YAxis.FontWeight = 'bold';


grid on;
box on;

if out2file
    outdir = sprintf('%s/%s/figures', params.io.output_dir, data.subject);
    if ~exist(outdir,'dir')
       mkdir(outdir); 
    end
    saveas(h, sprintf('%s/eeg_hilbert_envelopes.fig', outdir));
    saveas(h, sprintf('%s/eeg_hilbert_envelopes.png', outdir));
    close(h);
end


end