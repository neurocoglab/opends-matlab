function outcome = get_outcome( sim2track, interval ) 
% Determines the points outcome of this passing epoch

    idx1 = find(sim2track.matrix.TrackerTime > interval(1),1)-1;
    idx2 = find(sim2track.matrix.TrackerTime > interval(2),1);
    
    idx_r = find(sim2track.reward_times >= sim2track.matrix.TrackerTime(idx1) & ...
                 sim2track.reward_times <= sim2track.matrix.TrackerTime(idx2));

    if ~isempty(idx_r)
        % If more than one are in this interval, take the last value as
        % the outcome
        outcome = sim2track.reward_magnitudes(idx_r(end));
    else
        outcome = nan;
    end
             
end

