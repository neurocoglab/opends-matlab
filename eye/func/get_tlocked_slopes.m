function [slopes, timepoints] = get_tlocked_slopes( tlocked, params_event )
% Get slopes for each timelocked series, based on local min/max optima around the
% event

N_events = size(tlocked,1);
slopes = nan(N_events,1);
timepoints = nan(N_events,2);

idx_event = params_event.prepost(1);
dt = (params_event.prepost(2)+params_event.prepost(1)) / size(tlocked,2);

for i = 1 : N_events
    
    ts_i = tlocked(i,:);
    [~,idx_max] = findpeaks(ts_i, 'MinPeakDistance', params_event.slope_min_peak_dist, ...
                                   'MinPeakWidth', params_event.slope_min_peak_width);
    [~,idx_min] = findpeaks(-ts_i, 'MinPeakDistance', params_event.slope_min_peak_width, ...
                                   'MinPeakWidth', params_event.slope_min_peak_width);
    
    idx0 = find(idx_event-idx_min > 0 & idx_event-idx_min < params_event.slope_prepost(1), ...
                1, 'last');
    idx1 = find(idx_max-idx_event > 0 & idx_max-idx_event < params_event.slope_prepost(2), ...
                1, 'last');
    
    if ~isempty(idx0) && ~isempty(idx1)
        idx_min = idx_min(idx0);
        idx_max = idx_max(idx1);
        delta_t = (idx_max-idx_min) * dt;

        slopes(i) = (ts_i(idx_max) - ts_i(idx_min)) / delta_t;
        timepoints(i,:) = (-idx_event + [idx_min idx_max]) * dt;
    end
    
end




end