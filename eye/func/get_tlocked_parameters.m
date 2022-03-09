function tlock_params = get_tlocked_parameters( tlocked, params_event )
% Get slopes & amplitudes for each timelocked series, based on local min/max optima around the
% event

N_events = size(tlocked,1);
tlock_params = [];
tlock_params.max_dt = nan(N_events,1);
tlock_params.slopes = nan(N_events,1);
tlock_params.amplitudes = nan(N_events,1);
tlock_params.timepoints = nan(N_events,2);

idx_event = params_event.prepost(1);

% How far before event is a peak considered (ms)?
peak_window = 500;

for i = 1 : N_events
    
    if iscell(tlocked)
        ts_i = tlocked{i};
        ts_i = ts_i.avg;
    else
        ts_i = tlocked(i,:);
    end
    
    dt = (params_event.prepost(2)+params_event.prepost(1)) / length(ts_i);

    if params_event.tlock_params.smooth_span > 1
        ts_i = smooth(ts_i, params_event.tlock_params.smooth_span, 'lowess');
    end
    [~,idx_max] = findpeaks(ts_i, 'MinPeakDistance', params_event.tlock_params.min_peak_dist, ...
                                   'MinPeakWidth', params_event.tlock_params.min_peak_width);
    [~,idx_min] = findpeaks(-ts_i, 'MinPeakDistance', params_event.tlock_params.min_peak_width, ...
                                   'MinPeakWidth', params_event.tlock_params.min_peak_width);
    
    idx0 = find(idx_event-idx_min > 0 & idx_event-idx_min < params_event.tlock_params.prepost(1));
    idx1 = find(idx_max-idx_event > -peak_window & idx_max-idx_event < params_event.tlock_params.prepost(2));
    
    if isempty(idx0)
       idx_min = idx_event - params_event.tlock_params.prepost(1);
       idx0 = 1;
    end
   
    if ~isempty(idx1)
        % We have a valid peak, look for maximum temporal derivative in
        % preceding window
        t = ((1:length(ts_i))*dt - params_event.prepost(1));
        
        [ts_min,jj]= min(ts_i(idx_min(idx0)));
        idx_min = idx_min(idx0(jj));
        [ts_max,jj]= max(ts_i(idx_max(idx1)));
        idx_max = idx_max(idx1(jj));
        
        %if ts_max > ts_min
        delta_t = t(idx_max) - t(idx_min);
        delta_y = ts_i(idx_max) - ts_i(idx_min);
        
        d_ts = diff(ts_i(idx_min:idx_max));
        
        tlock_params.max_dt(i) = max(d_ts);
        tlock_params.slopes(i) = delta_y / delta_t;
        tlock_params.amplitudes(i) = delta_y;
        tlock_params.timepoints(i,:) = [t(idx_min),t(idx_max)];
        
    else
        wtf = 0;
        
    end
    
end




end