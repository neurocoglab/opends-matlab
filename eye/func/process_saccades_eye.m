function [ data ] = process_saccades_eye( data, params )
% PROCESS_SACCADES Searches eye tracking data (X and Y position) for
% saccades, defined by velocity thresholds specified in params

if isfield(data.eye, 'blinks')
    xx = data.eye.blinks.(params.eye.saccades.xpos_variable);
    yy = data.eye.blinks.(params.eye.saccades.ypos_variable);
else
    xx = data.eye.(params.eye.saccades.xpos_variable);
    yy = data.eye.(params.eye.saccades.ypos_variable);
end

data.eye.saccades = [];
data.eye.saccades.x_fixed = [];
data.eye.saccades.y_fixed = [];
data.eye.saccades.dx = [];
data.eye.saccades.dy = [];
data.eye.saccades.velocity = [];
data.eye.saccades.saccades = [];
data.eye.saccades.saccade_rate = [];

peaks = load(params.eye.saccades.peaks_file);
idxi = 1;

% tgaps = find(diff(data.eye.t) > 1/data.eye.Fs*1000);

for ii = 1 : size(data.eye.tgap,1)
%     idx1 = tgaps(ii); % data.tgap(ii,1);
    idx1 = data.eye.tgap(ii,1);
    
    xxi = xx(idxi:idx1);
    yyi = yy(idxi:idx1);
   
    data.eye.saccades = process_segment(idxi, xxi, yyi, peaks, data.eye.saccades);
    
    idxi = idx1 + 1;
    
end

if idxi < length(data.eye.t)
    xxi = xx(idxi:end);
    yyi = yy(idxi:end);
    data.eye.saccades = process_segment(idxi, xxi, yyi, peaks, data.eye.saccades);
elseif idxi == length(data.eye.t)
    % Degenerate case
    data.eye.saccades.x_fixed(end+1) = data.eye.saccades.x_fixed(end);
    data.eye.saccades.y_fixed(end+1) = data.eye.saccades.y_fixed(end);
    data.eye.saccades.saccade_rate(end+1) = data.eye.saccades.saccade_rate(end);
    data.eye.saccades.dx(end+1) = data.eye.saccades.dx(end);
    data.eye.saccades.dy(end+1) = data.eye.saccades.dy(end);
    data.eye.saccades.velocity(end+1) = data.eye.saccades.velocity(end);
end


    function [ sresults ] = process_segment(idxi, x, y, peaks, sresults)

        x = smooth(x,50,'sgolay');
        y = smooth(y,50,'sgolay');

        % Interpolate around large artifacts & eyeblinks
%         buffer = round(params.eye.saccades.artifact_buffer / data.Fs * 1000 / 2);

        % Convert to visual angle
        w_cm = params.eye.saccades.monitor_dims(1)/10;
        w_px = params.eye.saccades.monitor_res(1);
        d_cm = params.eye.saccades.distance/10;
        deg_per_px = radtodeg(atan2(0.5*w_cm, d_cm)) / (0.5*w_px);
        x = (x - w_px/2) * deg_per_px;

        h_cm = params.eye.saccades.monitor_dims(2)/10;
        h_px = params.eye.saccades.monitor_res(2);
        deg_per_px = radtodeg(atan2(0.5*h_cm, d_cm)) / (0.5*h_px);
        y = (y - h_px/2) * deg_per_px;

        sresults.x_fixed = [sresults.x_fixed;x];
        sresults.y_fixed = [sresults.y_fixed;y];

        mx = mean(x); vx = std(x);
        my = mean(y); vy = std(y);

        % Smooth position signal with smoothed boxcar shape
        for i = 1 : params.eye.saccades.nitr_smooth
            x = conv(x,peaks.skew,'same');
            y = conv(y,peaks.skew,'same');
        end

        x = (zscore(x) * vx) + mx;
        y = (zscore(y) * vy) + my;

        % Compute derivatives and velocity (vis. deg per sec)
        dx = [0;diff(x)] * data.eye.Fs;
        dy = [0;diff(y)] * data.eye.Fs;
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
                if velocity(i) > params.eye.saccades.velocity_thres
                    i_s = i;
                    issub = false;
                end
            else
                if velocity(i) < params.eye.saccades.velocity_thres
                    issub = true;
                    if i_s > 0
                        if i - i_s > params.eye.saccades.min_width && i - i_s < params.eye.saccades.max_width
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
        to_Hz = 1000 / params.eye.saccades.sr_window;
        Fs_ms = data.eye.Fs / 1000;
        hw = round(Fs_ms * params.eye.saccades.sr_window / 2);
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

