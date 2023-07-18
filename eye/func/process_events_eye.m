function [ results, summary ] = process_events_eye( params, data, results, summary )

if isempty(summary)
    summary = [];
    summary.overtake.tlocked = {};
    summary.overtake.tlocked_bl = {};
    summary.overtake.tlocked_bl2 = {};
    summary.overtake.diffs = {};
    summary.overtake.outcomes = {};

    summary.left_change.tlocked = {};
    summary.left_change.tlocked_bl = {};
    summary.left_change.tlocked_bl2 = {};
    summary.left_change.diffs = {};
    summary.left_change.outcomes = {};

    summary.right_change.tlocked = {};
    summary.right_change.tlocked_bl = {};
    summary.right_change.tlocked_bl2 = {};
    summary.right_change.diffs = {};
    summary.right_change.outcomes = {};
    
    % Process traffic decision button press events?
    if params.sim.events.traffic_decision.apply
        summary.traffic_decision.tlocked = {};
        summary.traffic_decision.tlocked_bl = {};
        summary.traffic_decision.tlocked_bl2 = {};
        summary.traffic_decision.confidence = {};
        summary.traffic_decision.correct = {};
        summary.traffic_decision.order = {};
        summary.traffic_decision.subjects = {};
    end
    
    % Process fixation onset events
    if params.sim.events.fixations.apply
        summary.fixations.aois = {};
        summary.fixations.tlocked = {};
        summary.fixations.tlocked_bl = {};
        summary.fixations.tlocked_bl2 = {};
        
    end

    summary.subjects = [];
    
end

if ~isfield(data.eye, 'blinks')
   warning( 'No eye blink correction has been performed; using raw data!' ); 
end

summary.subjects = [summary.subjects {data.subject}];

% Load sequence difficulty ratings
if params.eye.events.difficulty.apply
    seq_diff = readtable(sprintf('%s/%s/%s', params.io.input_dir, ...
                                             params.io.metadata_dir, ...
                                             params.sim.sequence_difficulty_file));
end

if isfield(data.eye, 'luminance') && ~data.eye.luminance.deficient
    idx_offset = params.eye.luminance.use_offset;
    if idx_offset < 1
       [~,idx_offset] = max(data.eye.luminance.r2); 
    end
    t_pd = data.eye.luminance.ts{idx_offset};
    pd = data.eye.luminance.diam{idx_offset};
else
    t_pd = data.eye.t;
    if isfield(data.eye, 'blinks')
        pd = data.eye.blinks.diam;
    else
        pd = data.eye.diam;
    end
end

if params.eye.events.zscore
    pdz = zscore(pd);
    pdz_raw = zscore(data.eye.diam);
else
    pdz = pd;
    pdz_raw = data.eye.diam;
end

if params.eye.events.smooth > 0
   pdz = smooth(pdz, params.eye.events.smooth);
   pdz_raw = smooth(pdz_raw, params.eye.events.smooth);
end

pdz = pdz(:);
t = data.eye.t;

% Get baselines
baseline = data.sim.sim2track.baseline;
idx_baseline = zeros(0,2);
t_baselines = [];
for i = 1 : size(baseline,1)
   
    % Map baseline times to ts indexes
    ti = baseline(i,1);
%     c = find(data.sim.sim2track.cycle_times > ti, 1);
%     if isempty(c); c = length(data.sim.sim2track.cycle_times)+1; end
    idx1 = find(t_pd < ti,1,'last');
    idx2 = idx1 + find(t_pd(idx1+1:end) > baseline(i,2),1,'first');
    
    if idx2 > idx1
        idx_baseline(end+1,:) = [idx1,idx2];
        t_baselines(end+1) = t_pd(idx1);
    end
    
end

N_rand = params.eye.events.random.N;

% Overtake events
events = data.sim.sim2track.overtake_times;
events = events(~isnan(events));

[tlocked,tstart] = get_tlocked(pdz, t_pd, events, params.eye.events.overtake);
% if params.eye.events.tlock_params.apply
%     [tlock_params] = get_tlocked_parameters(mean(tlocked), params.eye.events.overtake);
%     results.eye.events.overtake.tlock_params = tlock_params;
% end

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
results.eye.events.overtake.baseline_stats = baseline_stats;
results.eye.events.overtake.events = events;
results.eye.events.overtake.tlocked = tlocked;
results.eye.events.overtake.tlocked_bl = tlocked_bl;
results.eye.events.overtake.t = -params.eye.events.overtake.prepost(1):params.eye.events.overtake.prepost(2);
results.eye.events.overtake.t = results.eye.events.overtake.t/params.eye.Fs;
results.eye.events.overtake.diffs = zeros(length(events),1);
results.eye.events.overtake.outcomes = nan(length(events),1);

for i = 1 : length(events)
    if params.eye.events.difficulty.apply
        results.eye.events.overtake.diffs(i) = get_difficulty(seq_diff, data.sim.sim2track.matrix, ...
                                                              events(i), params.sim.lane_dist);                                             
    end
    if params.eye.events.outcomes.apply
        idxi = find(results.eye.epochs.overtake_intervals(:,1) <= events(i) & ...
                    results.eye.epochs.overtake_intervals(:,2) >= events(i));
        if ~isempty(idxi)
            results.eye.events.overtake.outcomes(i) = get_outcome(data.sim.sim2track, ...
                                                                  results.eye.epochs.overtake_intervals(idxi,:));
        end
    end
end

% Get corresponding random time-locked events for statistical comparison
N_event = length(events);
events = get_random_events( events, N_rand, [results.eye.events.overtake.t(1) results.eye.events.overtake.t(end)], ...
                            idx_baseline );
T = zeros(N_rand,N_event,length(results.eye.events.overtake.t));
T_params = [];
T_params.slope = zeros(N_rand,N_event);
T_params.amplitude = zeros(N_rand,N_event);
for i = 1 : N_rand
    tlocked = get_tlocked(pdz, t_pd, events(:,i), params.eye.events.overtake);
%     if params.eye.events.tlock_params.apply
%         tlock_params = get_tlocked_parameters(tlocked, params.eye.events.overtake);
%         T_params.slope(i,:) = tlock_params.slopes;
%         T_params.amplitude(i,:) = tlock_params.amplitudes;
%     end
    T(i,:,:) = tlocked;
end
results.eye.events.overtake.tlocked_bl2 = squeeze(nanmean(T,1));
% if params.eye.events.tlock_params.apply
%     results.eye.events.overtake.tlocked_params_bl2 = T_params;
% end

% Lane change left events
events = data.sim.sim2track.left_change_times;
events = events(~isnan(events));

[tlocked,tstart] = get_tlocked(pdz, t_pd, events, params.eye.events.left_change);
% if params.eye.events.tlock_params.apply
%     tlock_params = get_tlocked_parameters(tlocked, params.eye.events.left_change);
%     results.eye.events.left_change.tlock_params = tlock_params;
% end
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
results.eye.events.left_change.baseline_stats = baseline_stats;
results.eye.events.left_change.events = events;
results.eye.events.left_change.tlocked = tlocked;
results.eye.events.left_change.tlocked_bl = tlocked_bl;
results.eye.events.left_change.t = -params.eye.events.left_change.prepost(1):params.eye.events.left_change.prepost(2);
results.eye.events.left_change.t = results.eye.events.left_change.t/params.eye.Fs;
results.eye.events.left_change.diffs = zeros(length(events),1);
results.eye.events.left_change.outcomes = nan(length(events),1);

for i = 1 : length(events)
    if params.eye.events.difficulty.apply
        results.eye.events.left_change.diffs(i) = get_difficulty(seq_diff, data.sim.sim2track.matrix, ...
                                                                 events(i), params.sim.lane_dist);
    end
    if params.eye.events.outcomes.apply
        idxi = find(results.eye.epochs.overtake_intervals(:,1) <= events(i) & ...
                    results.eye.epochs.overtake_intervals(:,2) >= events(i));
        if ~isempty(idxi)
            results.eye.events.left_change.outcomes(i) = get_outcome(data.sim.sim2track, ...
                                                                     results.eye.epochs.overtake_intervals(idxi,:));
        end
    end
end

% Get corresponding random time-locked events for statistical comparison
N_event = length(events);
events = get_random_events( events, N_rand, [results.eye.events.left_change.t(1) results.eye.events.left_change.t(end)], ...
                            idx_baseline );
T = zeros(N_rand,N_event,length(results.eye.events.left_change.t));
T_params = [];
T_params.slope = zeros(N_rand,N_event);
T_params.amplitude = zeros(N_rand,N_event);
for i = 1 : N_rand
    tlocked = get_tlocked(pdz, t_pd, events(:,i), params.eye.events.left_change);
    T(i,:,:) = tlocked;
end
results.eye.events.left_change.tlocked_bl2 = squeeze(nanmean(T,1));
                                                 
% Lane change right events
events = data.sim.sim2track.right_change_times;
events = events(~isnan(events));

[tlocked,tstart] = get_tlocked(pdz, t_pd, events, params.eye.events.right_change);

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
results.eye.events.right_change.baseline_stats = baseline_stats;
results.eye.events.right_change.events = events;
results.eye.events.right_change.tlocked = tlocked;
results.eye.events.right_change.tlocked_bl = tlocked_bl;
results.eye.events.right_change.t = -params.eye.events.right_change.prepost(1):params.eye.events.right_change.prepost(2);
results.eye.events.right_change.t = results.eye.events.right_change.t/params.eye.Fs;
results.eye.events.right_change.diffs = zeros(length(events),1);
results.eye.events.right_change.outcomes = nan(length(events),1);

for i = 1 : length(events)
    if params.eye.events.difficulty.apply
        results.eye.events.right_change.diffs(i) = get_difficulty(seq_diff, data.sim.sim2track.matrix, ...
                                                                 events(i), params.sim.lane_dist);
    end
    if params.eye.events.outcomes.apply
        idxi = find(results.eye.epochs.overtake_intervals(:,1) <= events(i) & ...
                    results.eye.epochs.overtake_intervals(:,2) >= events(i));
        if ~isempty(idxi)
            results.eye.events.right_change.outcomes(i) = get_outcome(data.sim.sim2track, ...
                                                                      results.eye.epochs.overtake_intervals(idxi,:));
        end
    end
end

% Get corresponding random time-locked events for statistical comparison
N_event = length(events);
events = get_random_events( events, N_rand, [results.eye.events.right_change.t(1) results.eye.events.right_change.t(end)], ...
                            idx_baseline );
T = zeros(N_rand,N_event,length(results.eye.events.right_change.t));
T_params = [];
T_params.slope = zeros(N_rand,N_event);
T_params.amplitude = zeros(N_rand,N_event);
for i = 1 : N_rand
    tlocked = get_tlocked(pdz, t_pd, events(:,i), params.eye.events.right_change);
    T(i,:,:) = tlocked;
end
results.eye.events.right_change.tlocked_bl2 = squeeze(nanmean(T,1));


% Saccade offset events
svel = data.eye.saccades.saccades(:,4);
events = data.eye.saccades.saccades(svel>params.eye.events.saccades.vmin,3);
events=t(events);
tlocked = get_tlocked(pdz_raw, t, events', params.eye.events.saccades);
results.eye.events.saccades.events = events';
results.eye.events.saccades.tlocked = tlocked;
results.eye.events.saccades.t = -params.eye.events.saccades.prepost(1):params.eye.events.saccades.prepost(2);
results.eye.events.saccades.t = results.eye.events.saccades.t/params.eye.Fs;

% Blink events
%events = get_blink_offsets( preprocess.blink_ints );
events = get_blink_events( t, data.eye.blinks.blink_ints );
tlocked = get_tlocked(pdz_raw, t, events', params.eye.events.blinks);
results.eye.events.blinks.tlocked = tlocked;
results.eye.events.blinks.t = -params.eye.events.blinks.prepost(1):params.eye.events.blinks.prepost(2);
results.eye.events.blinks.t = results.eye.events.blinks.t/params.eye.Fs;

% Traffic decision events
if params.sim.events.traffic_decision.apply
    
    signals = {};
    signals.pdz = pdz;
    signals.t_pd = t_pd;
    signals.t_baselines = t_baselines;
    signals.idx_baseline = idx_baseline;
    
    results = process_traffic_decision_events( signals, data, results, params );
end


% Update summary
summary.overtake.tlocked = [summary.overtake.tlocked {results.eye.events.overtake.tlocked}];
summary.overtake.tlocked_bl = [summary.overtake.tlocked_bl {results.eye.events.overtake.tlocked_bl}];
summary.overtake.tlocked_bl2 = [summary.overtake.tlocked_bl2 {results.eye.events.overtake.tlocked_bl2}];
summary.overtake.diffs = [summary.overtake.diffs {results.eye.events.overtake.diffs}];
summary.overtake.outcomes = [summary.overtake.outcomes {results.eye.events.overtake.outcomes}];

summary.left_change.tlocked = [summary.left_change.tlocked {results.eye.events.left_change.tlocked}];
summary.left_change.tlocked_bl = [summary.left_change.tlocked_bl {results.eye.events.left_change.tlocked_bl}];
summary.left_change.tlocked_bl2 = [summary.left_change.tlocked_bl2 {results.eye.events.left_change.tlocked_bl2}];
summary.left_change.diffs = [summary.left_change.diffs {results.eye.events.left_change.diffs}];
summary.left_change.outcomes = [summary.left_change.outcomes {results.eye.events.left_change.outcomes}];

summary.right_change.tlocked = [summary.right_change.tlocked {results.eye.events.right_change.tlocked}];
summary.right_change.tlocked_bl = [summary.right_change.tlocked_bl {results.eye.events.right_change.tlocked_bl}];
summary.right_change.tlocked_bl2 = [summary.right_change.tlocked_bl2 {results.eye.events.right_change.tlocked_bl2}];
summary.right_change.diffs = [summary.right_change.diffs {results.eye.events.right_change.diffs}];
summary.right_change.outcomes = [summary.right_change.outcomes {results.eye.events.right_change.outcomes}];

% Traffic decisions
if params.sim.events.traffic_decision.apply
    if ~isempty(results.eye.events.traffic_decision)
        summary.traffic_decision.tlocked = [summary.traffic_decision.tlocked {results.eye.events.traffic_decision.tlocked}];
        summary.traffic_decision.tlocked_bl = [summary.traffic_decision.tlocked_bl {results.eye.events.traffic_decision.tlocked_bl}];
        summary.traffic_decision.tlocked_bl2 = [summary.traffic_decision.tlocked_bl2 {results.eye.events.traffic_decision.tlocked_bl2}];
        summary.traffic_decision.confidence = [summary.traffic_decision.confidence {results.eye.events.traffic_decision.confidence}];
        summary.traffic_decision.correct = [summary.traffic_decision.correct {results.eye.events.traffic_decision.correct}];
         summary.traffic_decision.order = [summary.traffic_decision.order {results.eye.events.traffic_decision.order}];
        if ~isfield(summary.traffic_decision, 't')
            summary.traffic_decision.t = results.eye.events.traffic_decision.t;
        end
        summary.traffic_decision.subjects = [summary.traffic_decision.subjects {data.subject}];
    end
end

if ~isfield(summary.overtake, 't')
    summary.overtake.t = results.eye.events.overtake.t;
    summary.left_change.t = results.eye.events.left_change.t;
    summary.right_change.t = results.eye.events.right_change.t;
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




