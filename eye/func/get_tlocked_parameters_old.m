function tlock_params = get_tlocked_parameters( tlocked, params_event )
% Get slopes & amplitudes for each timelocked series, based on local min/max optima around the
% event

N_events = size(tlocked,1);
tlock_params = [];
tlock_params.slopes = nan(N_events,1);
tlock_params.amplitudes = nan(N_events,1);
tlock_params.timepoints = nan(N_events,2);

idx_event = params_event.prepost(1);
dt = (params_event.prepost(2)+params_event.prepost(1)) / size(tlocked,2);

% How far before event is a peak considered (ms)?
peak_window = 500;

for i = 1 : N_events
    
    if iscell(tlocked)
        ts_i = tlocked{i};
        ts_i = ts_i.avg;
    else
        ts_i = tlocked(i,:);
    end

    if params_event.tlock_params.smooth_span > 1
        ts_i = smooth(ts_i, params_event.tlock_params.smooth_span, 'lowess');
    end
    [~,idx_max] = findpeaks(ts_i, 'MinPeakDistance', params_event.tlock_params.min_peak_dist, ...
                                   'MinPeakWidth', params_event.tlock_params.min_peak_width);
    [~,idx_min] = findpeaks(-ts_i, 'MinPeakDistance', params_event.tlock_params.min_peak_width, ...
                                   'MinPeakWidth', params_event.tlock_params.min_peak_width);
    
    idx0 = find(idx_event-idx_min > 0 & idx_event-idx_min < params_event.tlock_params.prepost(1));
    idx1 = find(idx_max-idx_event > -peak_window & idx_max-idx_event < params_event.tlock_params.prepost(2));
    
    if ~isempty(idx0) && ~isempty(idx1)
        [ts_min,jj]= min(ts_i(idx_min(idx0)));
        idx_min = idx_min(idx0(jj));
        [ts_max,jj]= max(ts_i(idx_max(idx1)));
        idx_max = idx_max(idx1(jj));
        
        %if ts_max > ts_min
            delta_t = (idx_max-idx_min) * dt;
            amplitude = ts_i(idx_max) - ts_i(idx_min);

            tlock_params.slopes(i) = amplitude / delta_t;
            tlock_params.amplitudes(i) = amplitude;
            tlock_params.timepoints(i,:) = (-idx_event + [idx_min idx_max]) * dt;
        %end

    end
    
end




end