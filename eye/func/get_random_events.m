function rand_events = get_random_events( events, N_rand, prepost, idx_bl )

    prepost(1) = -prepost(1);
    N_event = length(events);
    rand_events = nan(N_event, N_rand);

    for i = 1 : N_event
       
        idx_i = events(i);
        
        % Find closest baseline interval
        T = idx_bl(:,1)-idx_i;
        idx_cl1 = find(T>0,1);
        T = idx_i-idx_bl(:,2);
        idx_cl2 = find(T>0,1,'last');
        
        if isempty(idx_cl1) && isempty(idx_cl2)
           warning('No baseline found for events..');
           return;
        end
        
        if isempty(idx_cl1)
           idx_cl = idx_cl2;
        elseif isempty (idx_cl2)
           idx_cl = idx_cl1;
        elseif idx_bl(idx_cl1,1)-idx_i > idx_i-idx_bl(idx_cl2,2)
           idx_cl = idx_cl2;
        else
           idx_cl = idx_cl1;
        end
        
        % Assume baseline interval is longer than event window
        bl_int = idx_bl(idx_cl,:) - prepost;
        if any(bl_int<0)
            % ...
        else
            rand_events(i,:) = bl_int(1) + randi(bl_int(2),N_rand,1);
        end
    end
    
end