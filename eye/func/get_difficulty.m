function diff = get_difficulty( difficulty, matrix, time, lane_dist )
% Determines the difficulty assigned to this interval

    % Determine round + lane distance from matrix using time
    idx = find(matrix.TrackerTime > time,1);
    if isempty(idx)
       diff = nan;
       return 
    end
    if idx > 1, idx = idx - 1; end
    round = matrix.Cycle(idx); % {idx,3};
    dist = matrix.LaneDist(idx); % {idx,5};
    
    % Determine difficulty from seq_diff using round+dist
    idx = difficulty.Round == round;
    D = difficulty(idx,:);
    
    D_diff = D.Difficulty;
    D_start = D.From * lane_dist;
    D_end = D.To * lane_dist;
    
    idx = find(D_start <= dist & D_end >= dist);
    if isempty(idx)
       diff = nan;
       return 
    end
    
    % Use maximal difficulty for this interval
    diff = max(D_diff(idx));
    
    if diff < 3
        diff = 1;
    else
        diff = 2;
    end

end

