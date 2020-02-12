function [ results ] = process_saccades( results, data, params )
% PROCESS_SACCADES Searches eye tracking data (X and Y position) for
% saccades, defined by velocity thresholds specified in params

N = length(data.t);
xx = results.(params.xpos_variable);
yy = results.(params.ypos_variable);

results.saccades.t = data.t;
results.saccades.tgap = data.tgap;

results.saccades.x_fixed = [];
results.saccades.y_fixed = [];
results.saccades.dx = [];
results.saccades.dy = [];
results.saccades.velocity = [];
results.saccades.saccades = [];
results.saccades.saccade_rate = [];

peaks = load(params.peaks_file);
idxi = 1;

tgaps = find(diff(data.t) > 1/data.Fs*1000);

for ii = 1 : length(tgaps)
    idx1 = tgaps(ii); % data.tgap(ii,1);
    
    xxi = xx(idxi:idx1);
    yyi = yy(idxi:idx1);
   
    results.saccades = process_segment(idxi, xxi, yyi, peaks, results.saccades);
    
    idxi = idx1 + 1;
    
end

if idxi < length(data.t)
    xxi = xx(idxi:end);
    yyi = yy(idxi:end);
    results.saccades = process_segment(idxi, xxi, yyi, peaks, results.saccades);
elseif idxi == length(data.t)
    % Degenerate case
    results.saccades.x_fixed(end+1) = results.saccades.x_fixed(end);
    results.saccades.y_fixed(end+1) = results.saccades.y_fixed(end);
    results.saccades.saccade_rate(end+1) = results.saccades.saccade_rate(end);
    results.saccades.dx(end+1) = results.saccades.dx(end);
    results.saccades.dy(end+1) = results.saccades.dy(end);
    results.saccades.velocity(end+1) = results.saccades.velocity(end);
end


    function [ sresults ] = process_segment(idxi, x, y, peaks, sresults)

        x = smooth(x,50,'sgolay');
        y = smooth(y,50,'sgolay');

        % Interpolate around large artifacts & eyeblinks
%         buffer = round(params.artifact_buffer / data.Fs * 1000 / 2);

        % Convert to visual angle
        w_cm = params.monitor_dims(1)/10;
        w_px = params.monitor_res(1);
        d_cm = params.distance/10;
        deg_per_px = radtodeg(atan2(0.5*w_cm, d_cm)) / (0.5*w_px);
        x = (x - w_px/2) * deg_per_px;

        h_cm = params.monitor_dims(2)/10;
        h_px = params.monitor_res(2);
        deg_per_px = radtodeg(atan2(0.5*h_cm, d_cm)) / (0.5*h_px);
        y = (y - h_px/2) * deg_per_px;

        sresults.x_fixed = [sresults.x_fixed;x];
        sresults.y_fixed = [sresults.y_fixed;y];

        mx = mean(x); vx = std(x);
        my = mean(y); vy = std(y);

        % Smooth position signal with smoothed boxcar shape
        for i = 1 : params.nitr_smooth
            x = conv(x,peaks.skew,'same');
            y = conv(y,peaks.skew,'same');
        end

        x = (zscore(x) * vx) + mx;
        y = (zscore(y) * vy) + my;

        % Compute derivatives and velocity (vis. deg per sec)
        dx = [0;diff(x)] * data.Fs;
        dy = [0;diff(y)] * data.Fs;
        mdx = mean(dx); vdx = std(dx);
        mdy = mean(dy); vdy = std(dy);

        dx = conv(dx,peaks.pk,'same');
        dy = conv(dy,peaks.pk,'same');
        dx = conv(dx,peaks.pk3,'same');
        dy = conv(dy,peaks.pk3,'same');

        dx = (zscore(dx) * vdx) + mdx;
        dy = (zscore(dy) * vdy) + mdy;

        velocity = sqrt(dx.^2 + dy.^2) ;

        sresults.dx = [sresults.dx;dx];
        sresults.dy = [sresults.dy;dy];
        sresults.velocity = [sresults.velocity;velocity];

        % Threshold to identify saccades
        Nv = length(velocity);
        issub = true; i_s = 0;
        saccades = zeros(0,4);
        for i = 1 : Nv
            if issub    
                if velocity(i) > params.velocity_thres
                    i_s = i;
                    issub = false;
                end
            else
                if velocity(i) < params.velocity_thres
                    issub = true;
                    if i_s > 0
                        if i - i_s > params.min_width && i - i_s < params.max_width
                            [v_max,i_m] = max(velocity(i_s : i));
                            saccades(end+1,:) = [i_s, i, i_s + i_m, v_max];
                        end
                        i_s = 0;
                    end
                end
            end
        end

        saccades(:,1:3) = saccades(:,1:3) + idxi - 1;
        sresults.saccades = [sresults.saccades;saccades];

        % Saccade rate in Hz
        % Note: sr_window is specified in ms
        to_Hz = 1000 / params.sr_window;
        Fs_ms = data.Fs / 1000;
        hw = round(Fs_ms * params.sr_window / 2);
        saccade_rate = zeros(Nv,1);
        sacc = saccades(:,3);
        
        for i = 1 : Nv
           t0 = max(1, i-hw);
           t1 = min(Nv, i+hw);
           saccade_rate(i) = sum(sacc >= t0 & sacc <= t1);
        end

        saccade_rate = saccade_rate * to_Hz;
        saccade_rate = smooth(saccade_rate, 500);
        
        sresults.saccade_rate = [sresults.saccade_rate;saccade_rate];
    
    end

end

