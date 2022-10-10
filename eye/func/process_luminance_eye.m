function [ data ] = process_luminance_eye( data, params )
% Performs luminance correction on pupil diameter
%

    if isfield(data.eye, 'blinks')
        datain = data.eye.blinks;
    else
        datain = data.eye;
    end
    
    screen_dims = params.eye.saccades.monitor_res;
    w =  ceil(screen_dims(1) * params.eye.luminance.screen_pct/100);
    h = ceil(screen_dims(2) * params.eye.luminance.screen_pct/100);
    x_min = (screen_dims(1)-w)/2;
    x_max = x_min + w;
    y_min = (screen_dims(2)-h)/2;
    y_max = y_min + h;
    
    idx_center = datain.pos_x >= x_min & datain.pos_x <= x_max & ...
                 datain.pos_y >= y_min & datain.pos_y <= y_max;

    % Load luminance data
    lum_file = sprintf('%s/%s/%s/luminance.csv', params.io.input_dir, ...
                                                 params.eye.luminance.sub_dir, ...
                                                 data.subject);
    
    if ~exist(lum_file, 'file')
       return; 
    end
    
    lumtab = readtable(lum_file);
    y_lum = lumtab.luminance;
    idx_nan = ~isnan(y_lum);
    t_lum = lumtab.time(idx_nan);
    y_lum = y_lum(idx_nan);
    
    % Zero on first simulation event
    t_lum = t_lum - data.eye.t_start_sim;
    
    % Remove luminance outliers
    idx_keep = abs(zscore(y_lum))<params.eye.luminance.outlier_lim;
    y_lum = y_lum(idx_keep);
    t_lum = t_lum(idx_keep);
    ts_lum_orig = timeseries(y_lum, t_lum);
    
    clear lumtab t_lum y_lum;

    % Remove regressor outliers
    idx_keep = abs(zscore(datain.diam))<params.eye.luminance.outlier_lim;
    
    % Resample and match time series
    F = datain.diam(idx_keep & idx_center); % Fun fact: this transposes the vector
    ts_eye = timeseries(F(:), data.eye.t(idx_keep & idx_center));
    Fs = 1000/(ts_eye.Time(2)- ts_eye.Time(1));
    
    if params.eye.luminance.downsample > 0 && params.eye.luminance.downsample < Fs
        % Downsample to target frequency (prevents rank deficiency)
        step = round(Fs/params.eye.luminance.downsample);
        
        idx_ds = 1:step:length(ts_eye.Time);
        ts_eye = resample(ts_eye, ts_eye.Time(idx_ds));
        Fs2 = 1000/(ts_eye.Time(2)- ts_eye.Time(1));
        if params.general.debug
            fprintf(' [downsampled from %1.2f to %1.2f Hz]...', Fs, Fs2);
        end
    end
    
    ts_lum = resample(ts_lum_orig, ts_eye.Time);
    
    % Compute regression & residual error
    X = ts_lum.Data;
    y = ts_eye.Data;
    
    % Deal with NaNs
    idx_nan = isnan(X);
    if sum(idx_nan) > 0
        X = X(~idx_nan);
        y = y(~idx_nan);
    end
    if params.eye.luminance.smooth > 0
        X = smooth(X, params.eye.luminance.smooth);
        y = smooth(y, params.eye.luminance.smooth);
    end
    
    T_s = 1 / data.eye.Fs;
    N_offset = length(params.eye.luminance.offsets);
    data.eye.luminance = [];
    data.eye.luminance.lm = cell(N_offset,1);
    data.eye.luminance.diam = cell(N_offset,1);
    data.eye.luminance.rmse = zeros(N_offset,1);
    data.eye.luminance.r2 = zeros(N_offset,1);
    
    data.eye.luminance.deficient = false;
    
%     ts_eye = timeseries(datain.diam(:), double(data.eye.t_start) + data.eye.t);
    ts_eye = timeseries(datain.diam(:), data.eye.t);
    ts_lum = resample(ts_lum_orig, ts_eye.Time);
    idx_nan = find(isnan(ts_lum.Data));
    data.eye.luminance.idx = ~isnan(ts_lum.Data);
    if ~isempty(idx_nan)
        ts_eye = delsample(ts_eye,'Index',idx_nan);
        ts_lum = delsample(ts_lum,'Index',idx_nan);
    end
    
    data.eye.luminance.t = ts_lum.Time; % - double(data.eye.t_start);
    data.eye.luminance.y = ts_lum.Data;
    
    % Once for each time offset
    for i = 1 : N_offset
        
        offset = round(params.eye.luminance.offsets(i) / T_s);
        
        % Shift luminance by offset
        Xi = X(max(offset+1,1):min(end+offset,length(X)));
        yi = y(max(-offset+1,1):min(end-offset,length(X)));
        
        lm = fitlm(Xi, yi);
        if lm.Rsquared.Ordinary == 0
            data.eye.luminance.deficient = true;
        end
        
        data.eye.luminance.lm(i) = {lm};
        diam = datain.diam;
        idx_lum = max(-offset+1,1):min(length(ts_lum.Data)-offset,length(ts_lum.Data));
        ypred = lm.Coefficients.Estimate(1) + lm.Coefficients.Estimate(2) * ts_lum.Data(idx_lum);
        
        idx_eye = max(offset+1,1):min(length(ts_eye.Data)+offset,length(ts_eye.Data));
        resids = ts_eye.Data(idx_eye) - ypred;
        data.eye.luminance.rmse(i) = sqrt(mean(resids.^2));
        
        data.eye.luminance.r2(i) = corr(ypred,ts_eye.Data(idx_eye))^2;
        diam = mean(diam) + resids;
        data.eye.luminance.diam(i) = {diam};
        data.eye.luminance.ts(i) = {ts_eye.Time(idx_eye)}; %-double(data.eye.t_start)};
        data.eye.luminance.idx_offset(i) = {idx_eye};
        
    end

end

