function [ results ] = process_saccades( results, data, params )
% PROCESS_SACCADES Searches eye tracking data (X and Y position) for
% saccades, defined by velocity thresholds specified in params

% Convert to visual angle
h_cm = params.monitor_dims(2)/10;
h_px = params.monitor_res(2);
d_cm = params.distance/10;

deg_per_px = radtodeg(atan2(0.5*h_cm, d_cm)) / (0.5*h_px);

N = length(data.t);
x = results.(params.xpos_variable);
y = results.(params.ypos_variable);

x = x * deg_per_px;
y = y * deg_per_px;

% Remove tgaps
x = smooth(x,50,'sgolay');
y = smooth(y,50,'sgolay');

% Remove large artifacts
dx = zscore([0;diff(x)]);
dx=smooth(dx);
outliers = find(dx > params.artifact_thres);
dxdown = dx < -params.artifact_thres;

[vals,maxes,w,p] = findpeaks(dx,'MinPeakProminence',params.artifact_min_prominence);
maxes = maxes(vals > params.artifact_thres);

[vals,mins,w,p] = findpeaks(-dx,'MinPeakProminence',params.artifact_min_prominence);
mins = mins(vals > params.artifact_thres);

if ~isempty(maxes) && ~isempty(mins)
   
    for i = 1 : length(maxes)
       
       
    
    end
    
end

if ~isempty(outliers)
    do = [0;diff(outliers)];
   
    idx0 = 1;
    idx1 = 2;
    while idx1 < length(outliers)
        if do(idx1) > 1
            % Start of new outlier interval; find end point
            if idx1-idx0 > 10
                [~,idx2] = max(dx(idx0:idx1));
                idx2 = outliers(idx0 + idx2);
                idx3 = idx2+1;
                
                % Now find next minimum
                while idx3 <= N && idx3-idx2 < params.artifact_max_width && ~dxdown(idx3)
                    idx3 = idx3 + 1;
                end
                
                if idx3 <= N && ~dxdown(idx3)
                    % Have start and end indices; fix this artifact
                    idx2 = max(1,idx2-params.artifact_buffer);
                    idx3 = min(N,idx3+params.artifact_buffer);
                    x(idx2:idx3) = linterp([idx2 idx3],[x(idx2) x(idx3)],idx2:idx3);
                end
               
            end
            idx0 = idx1;
        end
        idx1 = idx1+1;
    end
end

% if ~isempty(artifacts)
%     for i = 1 : length(artifacts)
%         
%         xx = x(artifacts(i,1):artifacts(i,2));
%         
%         
%         
%         idx0 = outliers(i);
%         idx1 = outliers(i) + 1;
%         while idx1 < N && dx(idx1) > -params.artifact_thres
%             idx1 = idx1 + 1;
%         end
%         if idx1 < N && idx1 - idx0 < params.artifact_max_width
%            xs = x(max(1,idx0-1));
%            xe = x(idx1+1);
%            x(idx0:idx1) = linterp([idx0 idx1],[xs xe],idx0:idx1);
%         end
%     end
% end
% 
% dy = zscore([0;diff(y)]);
% outliers = find(dy > params.artifact_thres);
% if ~isempty(outliers)
%     for i = 1 : length(outliers)
%         idx0 = outliers(i);
%         idx1 = outliers(i) + 1;
%         while idx1 < N && dy(idx1) > -params.artifact_thres
%             idx1 = idx1 + 1;
%         end
%         if idx1 < N && idx1 - idx0 < params.artifact_max_width
%            ys = y(max(1,idx0-1));
%            ye = y(idx1+1);
%            y(idx0:idx1) = linterp([idx0 idx1],[ys ye],idx0:idx1);
%         end
%     end
% end

results.saccades.t = data.t;

% First pass - remove spikes
xma = smooth(x,50000);
xma = smooth(xma,1000);
yma = smooth(y,50000);
yma = smooth(yma,1000);

xx = abs(zscore(x));
yy = abs(zscore(y));
isgap = find(xx > params.spike_thres(1) | yy > params.spike_thres(1));
if ~isempty(isgap)
    x(isgap) = xma(isgap);
    y(isgap) = yma(isgap);
    isgap = isgap-1;
    if isgap(1) < 1, isgap = isgap(2:end); end
    x(isgap) = xma(isgap);
    y(isgap) = yma(isgap);
    isgap = isgap+2;
    if isgap(end) > length(x), isgap = isgap(1:end-1); end
    x(isgap) = xma(isgap);
    y(isgap) = yma(isgap);
end

% Second pass
xma = smooth(x,50000);
xma = smooth(xma,1000);
yma = smooth(y,50000);
yma = smooth(yma,1000);
results.saccades.x_ma = xma;
results.saccades.y_ma = yma;

%isgap = zeros(N,1);
x = x - xma;
x = x + mean(xma);
y = y - yma;
y = y + mean(yma);

% Remove spikes (set to mean)
% for i = 1 : 2
xx = abs(zscore(x));
yy = abs(zscore(y));
isgap = find(xx > params.spike_thres(2) | yy > params.spike_thres(2));
if ~isempty(isgap)
    x(isgap) = xma(isgap);
    y(isgap) = yma(isgap);
    isgap = isgap-1;
    if isgap(1) < 1, isgap = isgap(2:end); end
    x(isgap) = xma(isgap);
    y(isgap) = yma(isgap);
    isgap = isgap+2;
    if isgap(end) > length(x), isgap = isgap(1:end-1); end
    x(isgap) = xma(isgap);
    y(isgap) = yma(isgap);
end

% Third pass
xma = smooth(x,10000);
xma = smooth(xma,1000);
yma = smooth(y,10000);
yma = smooth(yma,1000);

x = x - xma;
x = x + mean(xma);
y = y - yma;
y = y + mean(yma);

xx = abs(zscore(x));
yy = abs(zscore(y));
isgap = find(xx > params.spike_thres | yy > params.spike_thres);
x(isgap) = xma(isgap);
y(isgap) = yma(isgap);
isgap = isgap-1;
if isgap(1) < 1, isgap = isgap(2:end); end
x(isgap) = xma(isgap);
y(isgap) = yma(isgap);
isgap = isgap+2;
if isgap(end) > length(x), isgap = isgap(1:end-1); end
x(isgap) = xma(isgap);
y(isgap) = yma(isgap);

% x = x(~isgap);
% y = y(~isgap);
% results.saccades.t = results.saccades.t(~isgap);
% 
% % Re-establish tgaps
dt = diff(results.saccades.t);
idx = dt > data.Fs + 1;
dt = dt(idx);
results.saccades.tgap = [find(idx), dt];

for i = 1 : length(results.saccades.tgap)
   j = results.saccades.tgap(i,1);
   radius = 50;
   j0 = max(j - radius,1);
   j1 = min(j + radius,length(x));
   x(j0:j) = x(j0);
   x(j+1:j1) = x(j1);
   y(j0:j) = y(j0);
   y(j+1:j1) = y(j1);
end

for i = 1 : 2
    x = smooth(x,150,'loess');
    y = smooth(y,150,'loess');
end

% For X and Y position time series, get first derivative (velocity)
dx = [0;diff(x)];
dy = [0;diff(y)];

% isoutlier = abs(dx) > params.diff_thres | abs(dy) > params.diff_thres;
% isoutlier(results.saccades.tgap) = true;
% isoutlier(results.saccades.tgap+1) = true;
% dx(isoutlier) = 0;
% dy(isoutlier) = 0;

for i = 1 : length(results.saccades.tgap)
   j = results.saccades.tgap(i,1);
   gap = results.saccades.tgap(i,2);
   if gap > 50
       j0 = max(1,j-30); j1 = min(N,j+30); 
       dx(j0:j1) = 0;
       dy(j0:j1) = 0;
   end
end

velocity = sqrt(dx.^2+dy.^2);

results.saccades.dx = dx;
results.saccades.dy = dy;
results.saccades.velocity = velocity;

% Where velocity > threshold, define saccade
% Intervals between saccades are defined as fixation periods
results.saccades.x_fixed = x;
results.saccades.y_fixed = y;

% Identify saccades
N = length(velocity);
issub = true; i_s = 0;
saccades = zeros(0,4);
for i = 1 : N
    if issub    
        if velocity(i) > params.velocity_thres
            i_s = i;
            issub = false;
        end
    else
        if velocity(i) < params.velocity_thres
            issub = true;
            if i_s > 0
                [v_max,i_m] = max(velocity(i_s : i));
                if i_m - i_s > params.min_width
                    saccades(end+1,:) = [i_s, i, i_s + i_m, v_max];
                end
                i_s = 0;
            end
        end
    end
end

results.saccades.saccades = saccades;

% Saccade rate in Hz
hw = params.sr_window / 2;
saccade_rate = zeros(N,1);
sacc = saccades(:,3);
t = results.saccades.t;

w1 = 1; w2 = 1;

for i = 1 : N
   while w1 < N && t(w1) < t(i) - hw
       w1 = w1 + 1;
   end
   while w2 < N && t(w2) < t(i) + hw
       w2 = w2 + 1;
   end
   
    interval = (sacc >= w1 & sacc <= w2);
    saccade_rate(i) = sum(interval) / (t(w2) - t(w1));
   
end

saccade_rate = saccade_rate * 1000;
saccade_rate = smooth(saccade_rate, 500);
results.saccades.saccade_rate = saccade_rate;

end

