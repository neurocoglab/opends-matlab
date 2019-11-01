function [] = plot_eeg_blinks(eye_data, eeg_data, channels, stdev, to_file)

% Plot EEG and blinks together
% Time variables must have already been aligned

if nargin <5
    to_file = [];
end

idx_channels = [];
if nargin > 2
    for i = 1 : length(channels)
        idx_channels = [idx_channels find(strcmp(eeg_data.label,channels{i}))];
    end
end

if nargin < 4
   stdev = 3; 
end

if ~isempty(to_file)
    h = figure('visible','off');
else
    h = figure;
end
h.Color = [1 1 1];
x_eye = seconds(eye_data.t / 1000);

t_end = eye_data.sim2track.simended_time / 1000;
idx_keep = find(x_eye <= seconds(t_end));
x_eye = x_eye(idx_keep);

N_e = length(eye_data.blinks);

idx_sacc = eye_data.saccades.saccades(:,1:2);
idx_sacc = idx_sacc(~any(idx_sacc>max(idx_keep),2),:);
for i = 1 : size(idx_sacc,1)
    xs = x_eye(idx_sacc(i,:));
    hh = fill([xs(1) xs(1) xs(2) xs(2)],[-1000 1000 1000 -1000],[0.5 0.1 0.1]);     
    hh.EdgeColor = 'none';
    alpha(hh, 0.3);
    hold on;
end

% idx_blink = find(eye_data.blinks);
% y_blink = ones(length(idx_blink),1) * 1000;
% x_blink = x_eye(idx_blink);
% hh = stem(x_blink, y_blink, 'k');
% hh.LineWidth=2;
% hh = stem(x_blink, -y_blink, 'k');
% hh.LineWidth=2;

for i = 1 : length(eye_data.blink_ints)
    ints_i = eye_data.blink_ints{i};
    ints_i = ints_i(~any(ints_i>max(idx_keep),2),:);
    for j = 1 : size(ints_i,1)
        xs = [x_eye(ints_i(j,1)) x_eye(ints_i(j,1)+ints_i(j,2))];
        hh = fill([xs(1) xs(1) xs(2) xs(2)],[-1000 1000 1000 -1000],[0.1 0.1 0.5]);     
        hh.EdgeColor = 'none';
        alpha(hh, 0.3);
        hold on;
    end
end
    
L = 0;
lens = zeros(length(eeg_data.time),1);
for c = 1 : length(eeg_data.time)
    lens(c) = length(eeg_data.time{c});
    L = L + lens(c) + 1;
end
x_eeg = duration(nan(L,3));
y_eeg = zeros(length(idx_channels),L);
idx = 1;
for c = 1 : length(eeg_data.time)
    x_eeg(idx:idx+lens(c)-1) = seconds(eeg_data.time{c});
    x_eeg(idx+lens(c)) = x_eeg(idx+lens(c)-1) + seconds(0.00001);
    yy = eeg_data.trial{c};
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
gain = itr / stdev; % Std devs per lane

for i = 1 : N_ch
   y_eeg(i,:) = zscore_nan(y_eeg(i,:)) * gain - (i-1)*itr;
end

hh = plot(x_eeg, y_eeg);
ypos = zeros(N_ch,1);
ystr = cell(N_ch,1);
for i = 1 : N_ch
   ypos(i) = nanmean(y_eeg(i,:));
   ystr(i) = eeg_data.label(idx_channels(i));
end
ax = gca;

y_posx = 150 + zscore(eye_data.diam_left(idx_keep)) * 10;
hh = plot(x_eye, y_posx, 'b', 'LineWidth', 2);

y_posx = 110 + zscore(eye_data.pos_left_x(idx_keep)) * 10;
hh = plot(x_eye, y_posx, 'Color', [0 .3 0], 'LineWidth', 2);

y_posy = 75 + zscore(eye_data.pos_left_y(idx_keep)) * 10;
hh = plot(x_eye, y_posy, 'm', 'LineWidth', 2);

xlim(seconds([0 10]));
ylim([-300 200]);

xlabel('Time (mm:ss)');
ax.XAxis.Label.FontSize = 20;

ax.YTick=flip(ypos);
ax.YTickLabel=flip(ystr);
ax.XAxis.TickLabelFormat = 'mm:ss';

ax.FontSize=17;

if ~isempty(to_file)
    saveas(h,to_file);
    close(h);
end

    
function Z = zscore_nan(X)
        
    Z = (X-nanmean(X)) ./ nanstd(X);

end



end