function [ results ] = process_luminance( results, data, params, out2file )
% Performs luminance correction on pupil diameter
%
    results.luminance.deficient = false;

    if nargin < 4
        out2file = false;
    end
    
    screen_dims = params.saccades.monitor_res;
    w =  ceil(screen_dims(1) * params.luminance.screen_pct/100);
    h = ceil(screen_dims(2) * params.luminance.screen_pct/100);
    x_min = (screen_dims(1)-w)/2;
    x_max = x_min + w;
    y_min = (screen_dims(2)-h)/2;
    y_max = y_min + h;
    
    idx_center = results.pos_left_x >= x_min & results.pos_left_x <= x_max & ...
                 results.pos_left_y >= y_min & results.pos_left_y <= y_max;

    outdir = sprintf('%s/%s/%s/figures', params.root_dir, params.output_dir, data.subject);
    if ~exist(outdir,'dir')
       mkdir(outdir); 
    end

    % Load luminance data
    lum_file = sprintf('%s/%s/%s/luminance.csv', params.root_dir, params.luminance.dir, data.subject);
    lumtab = readtable(lum_file);
    y_lum = lumtab.luminance;
    idx_nan = ~isnan(y_lum);
    t_lum = lumtab.time(idx_nan);
    y_lum = y_lum(idx_nan);
    
    % Remove outliers
    idx_keep = abs(zscore(y_lum))<params.luminance.outlier_lim;
    y_lum = y_lum(idx_keep);
    t_lum = t_lum(idx_keep);
    ts_lum_orig = timeseries(y_lum, t_lum);
    
    clear lumtab t_lum y_lum;

    % Remove outliers first
    idx_keep = abs(zscore(results.diam_left))<params.luminance.outlier_lim;
    
    % Resample and match time series
    ts_eye = timeseries(results.diam_left(idx_keep & idx_center), ...
                        double(data.eye.t_start) + results.t(idx_keep & idx_center));
    if params.luminance.downsample > 1
        % Downsample by factor (prevents rank deficiency)
        idx_ds = 1:params.luminance.downsample:length(ts_eye.Time);
        ts_eye = resample(ts_eye, ts_eye.Time(idx_ds));
    end
    
    ts_lum = resample(ts_lum_orig, ts_eye.Time);
    
    % Compute regression & residual error
    X = ts_lum.Data; %(idx_keep & idx_center);
    y = ts_eye.Data; %(idx_keep & idx_center);
    
    % Deal with NaNs
    idx_nan = isnan(X);
    if sum(idx_nan) > 0
        X = X(~idx_nan);
        y = y(~idx_nan);
    end
    if params.luminance.smooth > 0
        X = smooth(X, params.luminance.smooth);
        y = smooth(y, params.luminance.smooth);
    end
    
    T_s = 1 / data.eye.Fs;
    N_offset = length(params.luminance.offsets);
    results.luminance.lm = cell(N_offset,1);
    results.luminance.diam_left = cell(N_offset,1);
    results.luminance.rmse = zeros(N_offset,1);
    results.luminance.r2 = zeros(N_offset,1);
    
    results.deficient = false;
    
    ts_eye = timeseries(results.diam_left, double(data.eye.t_start) + results.t);
    ts_lum = resample(ts_lum_orig, ts_eye.Time);
    idx_nan = find(isnan(ts_lum.Data));
    results.luminance.idx = ~isnan(ts_lum.Data);
    if ~isempty(idx_nan)
        ts_eye = delsample(ts_eye,'Index',idx_nan);
        ts_lum = delsample(ts_lum,'Index',idx_nan);
    end
    
    results.luminance.t = ts_lum.Time - double(results.t_start);
    results.luminance.y = ts_lum.Data;
    
    % Once for each time offset
    for i = 1 : N_offset
        
        offset = round(params.luminance.offsets(i) / T_s);
        
        % Shift luminance by offset
        Xi = X(max(offset+1,1):min(end+offset,length(X)));
        yi = y(max(-offset+1,1):min(end-offset,length(X)));
        
        lm = fitlm(Xi, yi);
        if lm.Rsquared.Ordinary == 0
            results.luminance.deficient = true;
        end
        
%         coefs = [ones(length(Xi),1),Xi]\yi;
        results.luminance.lm(i) = {lm};
        diam_left = results.diam_left;
        idx_lum = max(-offset+1,1):min(length(ts_lum.Data)-offset,length(ts_lum.Data));
        ypred = lm.Coefficients.Estimate(1) + lm.Coefficients.Estimate(2) * ts_lum.Data(idx_lum);
%         ypred = coefs(1) + coefs(2) * ts_lum.Data(idx_keep);
        
        idx_eye = max(offset+1,1):min(length(ts_eye.Data)+offset,length(ts_eye.Data));
        resids = ts_eye.Data(idx_eye) - ypred;
        results.luminance.rmse(i) = sqrt(mean(resids.^2));
        
        results.luminance.r2(i) = corr(ypred,ts_eye.Data(idx_eye))^2;
        diam_left = mean(diam_left) + resids;
        results.luminance.diam_left(i) = {diam_left};
        results.luminance.ts(i) = {ts_eye.Time(idx_eye)-double(results.t_start)};
        results.luminance.idx_offset(i) = {idx_eye};
        
    end

end

