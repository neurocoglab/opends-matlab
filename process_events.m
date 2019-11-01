function [ results ] = process_events( preprocess, data, results, params, preproc )

% Load sequence difficulty ratings
seq_diff = load(params.sequence_difficulty_file);

if isfield(preprocess, 'luminance') && ~preprocess.luminance.deficient
    idx_offset = preproc.params.luminance.use_offset;
    if idx_offset < 1
       [~,idx_offset] = max(preprocess.luminance.r2); 
    end
    t_pd = preprocess.luminance.ts{idx_offset};
    pd = preprocess.luminance.diam_left{idx_offset};
else
    t_pd = data.eye.t;
    pd = data.eye.diam_left;
end

if params.events.zscore
    pdz = zscore(smooth(pd,500));
    pdz_raw = zscore(smooth(data.eye.diam_left,500));
else
    pdz = smooth(pd,500);
    pdz_raw = smooth(data.eye.diam_left,500);
end

t = data.eye.t;

% Get baselines
baseline = preprocess.sim2track.baseline;
idx_baseline = zeros(0,2);
t_baselines = [];
for i = 1 : size(baseline,1)
   
    % Map baseline times to ts indexes
    ti = baseline(i,1);
    c = find(preprocess.sim2track.cycle_times > ti, 1);
    if isempty(c); c = length(preprocess.sim2track.cycle_times)+1; end
    idx1 = find(t_pd < ti,1,'last');
    idx2 = idx1 + find(t_pd(idx1+1:end) > baseline(i,2),1,'first');
    
    if idx2 > idx1
        idx_baseline(end+1,:) = [idx1,idx2];
        t_baselines(end+1) = t_pd(idx1);
    end
    
end


% Overtake outcomes


% Overtake events
events = preprocess.sim2track.overtake_times;

[tlocked,tstart] = get_tlocked(pdz, t_pd, events, params.events.overtake);
clear baseline_stats;
baseline_stats.mean = nan(length(tstart),1);
baseline_stats.std = nan(length(tstart),1);
tlocked_bl = zeros(size(tlocked));
for i = 1 : length(tstart)
   idx_bl = find(t_baselines < tstart(i),1,'last');
   if ~isempty(idx_bl)
      xx = pdz(idx_baseline(idx_bl,1):idx_baseline(idx_bl,2));
      baseline_stats.mean(i) = mean(xx);
      baseline_stats.std(i) = std(xx);
      tlocked_bl(i,:) = tlocked(i,:) - baseline_stats.mean(i);
   end
end
results.overtake.baseline_stats = baseline_stats;
results.overtake.events = events;
results.overtake.tlocked = tlocked;
results.overtake.tlocked_bl = tlocked_bl;
results.overtake.t = -params.events.overtake.prepost(1):params.events.overtake.prepost(2);
results.overtake.t = results.overtake.t/preproc.params.Fs;
results.overtake.diffs = zeros(length(events),1);
results.overtake.outcomes = nan(length(events),1);

for i = 1 : length(events)
    results.overtake.diffs(i) = get_difficulty(seq_diff.difficulty, preprocess.sim2track.matrix, ...
                                               events(i), params.lane_dist);
    idxi = find(results.epochs.overtake_intervals(:,1) <= events(i) & ...
                results.epochs.overtake_intervals(:,2) >= events(i));
    if ~isempty(idxi)
        results.overtake.outcomes(i) = get_outcome(preprocess.sim2track, ...
                                                   results.epochs.overtake_intervals(idxi,:), ...
                                                   params.lane_dist);
    end
end

% Lane change left events
events = preprocess.sim2track.left_change_times;

[tlocked,tstart] = get_tlocked(pdz, t_pd, events, params.events.left_change);
clear baseline_stats;
baseline_stats.mean = nan(length(tstart),1);
baseline_stats.std = nan(length(tstart),1);
tlocked_bl = zeros(size(tlocked));
for i = 1 : length(tstart)
   idx_bl = find(t_baselines < tstart(i),1,'last');
   if ~isempty(idx_bl)
      xx = pdz(idx_baseline(idx_bl,1):idx_baseline(idx_bl,2));
      baseline_stats.mean(i) = mean(xx);
      baseline_stats.std(i) = std(xx);
      tlocked_bl(i,:) = tlocked(i,:) - baseline_stats.mean(i);
   end
end
results.left_change.baseline_stats = baseline_stats;
results.left_change.events = events;
results.left_change.tlocked = tlocked;
results.left_change.tlocked_bl = tlocked_bl;
results.left_change.t = -params.events.left_change.prepost(1):params.events.left_change.prepost(2);
results.left_change.t = results.left_change.t/preproc.params.Fs;
results.left_change.diffs = zeros(length(events),1);
results.left_change.outcomes = nan(length(events),1);

for i = 1 : length(events)
    results.left_change.diffs(i) = get_difficulty(seq_diff.difficulty, preprocess.sim2track.matrix, ...
                                                  events(i), params.lane_dist);
    idxi = find(results.epochs.overtake_intervals(:,1) <= events(i) & ...
                results.epochs.overtake_intervals(:,2) >= events(i));
    if ~isempty(idxi)
        results.left_change.outcomes(i) = get_outcome(preprocess.sim2track, ...
                                                   results.epochs.overtake_intervals(idxi,:), ...
                                                   params.lane_dist);
    end
end

% Lane change right events
events = preprocess.sim2track.right_change_times;

[tlocked,tstart] = get_tlocked(pdz, t_pd, events, params.events.right_change);
clear baseline_stats;
baseline_stats.mean = nan(length(tstart),1);
baseline_stats.std = nan(length(tstart),1);
tlocked_bl = zeros(size(tlocked));
for i = 1 : length(tstart)
   idx_bl = find(t_baselines < tstart(i),1,'last');
   if ~isempty(idx_bl)
      xx = pdz(idx_baseline(idx_bl,1):idx_baseline(idx_bl,2));
      baseline_stats.mean(i) = mean(xx);
      baseline_stats.std(i) = std(xx);
      tlocked_bl(i,:) = tlocked(i,:) - baseline_stats.mean(i);
   end
end
results.right_change.baseline_stats = baseline_stats;
results.right_change.events = events;
results.right_change.tlocked = tlocked;
results.right_change.tlocked_bl = tlocked_bl;
results.right_change.t = -params.events.right_change.prepost(1):params.events.right_change.prepost(2);
results.right_change.t = results.right_change.t/preproc.params.Fs;
results.right_change.diffs = zeros(length(events),1);
results.right_change.outcomes = nan(length(events),1);

for i = 1 : length(events)
    results.right_change.diffs(i) = get_difficulty(seq_diff.difficulty, preprocess.sim2track.matrix, ...
                                                  events(i), params.lane_dist);
    idxi = find(results.epochs.overtake_intervals(:,1) <= events(i) & ...
                results.epochs.overtake_intervals(:,2) >= events(i));
    if ~isempty(idxi)
        results.right_change.outcomes(i) = get_outcome(preprocess.sim2track, ...
                                                   results.epochs.overtake_intervals(idxi,:), ...
                                                   params.lane_dist);
    end
end

% Saccade offset events
svel = preprocess.saccades.saccades(:,4);
events = preprocess.saccades.saccades(svel>params.events.saccades.vmin,3);
events=t(events);
tlocked = get_tlocked(pdz_raw, t, events', params.events.saccades);
results.saccades.events = events';
results.saccades.tlocked = tlocked;
results.saccades.t = -params.events.saccades.prepost(1):params.events.saccades.prepost(2);
results.saccades.t = results.saccades.t/preproc.params.Fs;

% Blink events
%events = get_blink_offsets( preprocess.blink_ints );
events = get_blink_events( t, preprocess.blink_ints );
tlocked = get_tlocked(pdz_raw, t, events', params.events.blinks);
results.blinks.tlocked = tlocked;
results.blinks.t = -params.events.blinks.prepost(1):params.events.blinks.prepost(2);
results.blinks.t = results.blinks.t/preproc.params.Fs;

% Permutations
% if params.events.perms.N > 0
% 
%     results.perms.t = -params.events.perms.prepost(1):params.events.perms.prepost(2);
%     results.perms.t = results.perms.t/preproc.params.Fs;
% 
%     for p = 1 : params.events.perms.N
%         events = t(randsample(length(t),params.events.perms.M));
%         tlocked_p = get_tlocked(pdz, t, events', params.events.perms);
%         if p == 1
%             nans_i = isnan(tlocked_p);
%             tlocked_p(nans_i) = 0;
%             tlocked = tlocked_p;
%             nan_vals = tlocked * 0;
%             nan_vals(nans_i) = 1;
%         else
%             nans_i = isnan(tlocked_p);
%             tlocked_p(nans_i) = 0;
%             tlocked = tlocked + tlocked_p;
%             nan_vals(nans_i) = nan_vals(nans_i) + 1;
%         end
%     end
% 
%     nan_vals = params.events.perms.N - nans_i;
%     results.perms.tlocked = tlocked ./ nan_vals;
% 
%     for p = 1 : params.events.perms.N
%         events = t(randsample(length(t),params.events.perms.M));
%         tlocked_p = get_tlocked(pdz_raw, t, events', params.events.perms);
%         if p == 1
%             nans_i = isnan(tlocked_p);
%             tlocked_p(nans_i) = 0;
%             tlocked = tlocked_p;
%             nan_vals = tlocked * 0;
%             nan_vals(nans_i) = 1;
%         else
%             nans_i = isnan(tlocked_p);
%             tlocked_p(nans_i) = 0;
%             tlocked = tlocked + tlocked_p;
%             nan_vals(nans_i) = nan_vals(nans_i) + 1;
%         end
%     end
% 
%     nan_vals = params.events.perms.N - nans_i;
%     results.perms.tlocked_raw = tlocked ./ nan_vals;
% 
% end

end

function [ tlocked, tstart ] = get_tlocked(pdz, t, events, params)

    tlocked = zeros(length(events),sum(params.prepost)+1);
    tstart = zeros(length(events),1);
    N=length(t);
        
    for i = 1 : length(events)

        idx = find(t>events(i),1,'first');
        if isempty(idx)
            warning('Event %d (%d) is beyond time vector.', i, events(i));
            break;
        end
        pdi = pdz(idx);
        pre = idx - params.prepost(1);
        if pre < 1
            pad0 = 1-pre;
            pre = 1;
        else
            pad0 = 0;
        end
        post = idx + params.prepost(2);
        if post > N
            pad1 = post-N;
            post = N;
        else
            pad1 = 0;
        end
        
        xx = pdz(pre:post);

        T = [nan(pad0,1);xx;nan(pad1,1)];
        tlocked(i,:) = T;
        tstart(i) = t(pre+pad0);

    end

end

% Blink is the center of an interval
function blink_events = get_blink_events( t, blink_ints )

    blink_events = [];

    for i = 1 : length(blink_ints)
        ints = blink_ints{i};
        for j = 1 : size(ints,1)
            idx = ints(j,1) + round(ints(j,2) / 2);
            blink_events = [blink_events t(idx)];
        end
    end

end

function offsets = get_blink_offsets( blink_ints )

    offsets = zeros(0,2);
    
    for i = 1 : length(blink_ints)
       offsets = [offsets; blink_ints{i}]; 
    end
    
    offsets = sum(offsets,2);

end

% Determines the difficulty assigned to this interval
function diff = get_difficulty( difficulty, matrix, time, lane_dist )

    % Determine round + lane distance from matrix using time
    idx = find(cell2mat(matrix(:,2))>time,1);
    if isempty(idx)
       diff = nan;
       return 
    end
    if idx > 1, idx = idx - 1; end
    round = matrix{idx,3};
    dist = matrix{idx,5};
    
    % Determine difficulty from seq_diff using round+dist
    idx = cell2mat(difficulty(:,1)) == round;
    D = difficulty(idx,:);
    
    D_diff = cell2mat(D(:,6));
    D_start = cell2mat(D(:,4)) * lane_dist;
    D_end = cell2mat(D(:,5)) * lane_dist;
    
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

% Determines the points outcome of this passing epoch
function outcome = get_outcome( sim2track, interval, lane_dist ) 

    idx1 = find(cell2mat(sim2track.matrix(:,2))>interval(1),1)-1;
    idx2 = find(cell2mat(sim2track.matrix(:,2))>interval(2),1);
    
    idx_r = find(sim2track.reward_times >= sim2track.matrix{idx1,2} & ...
                 sim2track.reward_times <= sim2track.matrix{idx2,2});

    if ~isempty(idx_r)
        % If more than one are in this interval, take the last value as
        % the outcome
        outcome = sim2track.reward_magnitudes(idx_r(end));
    else
        outcome = nan;
    end
             
end