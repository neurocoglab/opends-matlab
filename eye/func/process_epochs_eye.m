function [ results ] = process_epochs( preprocess, data, params, preproc )
% Summarize data over epochs
% Epochs:   baseline, non-baseline, passing events
% Data:     pupil diameter, saccade rate, blink rate

% Load sequence difficulty ratings
seq_diff = load(params.sequence_difficulty_file);
% D = cell2mat(seq_diff.difficulty(:,find(strcmp(seq_diff.hdr,'Difficulty'))));
diff_levels = [1 2]; % unique(D);

% Baseline
baseline = preprocess.sim2track.baseline;
 
if isfield(preprocess, 'luminance') && ~preprocess.luminance.deficient
    idx_offset = preproc.params.luminance.use_offset;
    if idx_offset < 1
       [~,idx_offset] = max(preprocess.luminance.r2); 
    end
    t_eye = preprocess.luminance.ts{idx_offset};
    pd = preprocess.luminance.diam_left{idx_offset};
else
    t_eye = data.eye.t;
    pd = data.eye.diam_left;
end

N = length(t_eye);
idx_baseline = false(N,1);
baseline_intervals = zeros(0,2);
baseline_cycles = zeros(N,1);

for i = 1 : size(baseline,1)
   
    % Map baseline times to ts indexes
    ti = baseline(i,1);
    c = find(preprocess.sim2track.cycle_times > ti, 1);
    if isempty(c); c = length(preprocess.sim2track.cycle_times)+1; end
    idx1 = find(t_eye < ti,1,'last');
    idx2 = idx1 + find(t_eye(idx1+1:end) > baseline(i,2),1,'first');
    
    if idx2 > idx1
        idx_baseline(idx1:idx2) = true;
        baseline_cycles(idx1:idx2) = c;
        try
        baseline_intervals(end+1,:) = [idx1 idx2];
        catch
           a=0; 
        end
    end
    
end

results.idx_baseline = idx_baseline;

% Passing
passing_intervals = zeros(0,2);
passing_times = zeros(0,2);
idx_passing = false(N,1);
passing_diff = zeros(N,1);
passing_outcomes = zeros(N,1);
passing_cycles = zeros(N,1);

max_interval = 30000; % Half a minute
        
left = [preprocess.sim2track.left_change_times, ...
        true(length(preprocess.sim2track.left_change_times),1)];
right = [preprocess.sim2track.right_change_times, ...
         false(length(preprocess.sim2track.right_change_times),1)];

left_right = [left;right];
[~,idx] = sort(left_right(:,1));
left_right = left_right(idx,:);

results.overtake_intervals = zeros(0,2);
results.overtake_outcomes = [];

this_left = -1;
this_right = -1;
window = round(params.epochs.overtake_window * data.eye.Fs);
for j = 1 : length(left_right)
   if left_right(j,2)
      % Is change to left 
      this_left = left_right(j,1);
   else
      % Is change to right
      this_right = left_right(j,1);
      if this_left > 0 && this_right-this_left < max_interval
          % Valid passing segment, add
          results.overtake_intervals(end+1,:) = [this_left this_right];
          
          c = find(preprocess.sim2track.cycle_times > this_left, 1);
          if isempty(c); c = length(preprocess.sim2track.cycle_times)+1; end
          
          % Map lane change times to ts indexes
          idx1 = find(t_eye < this_left,1,'last');
          ii = find(t_eye(idx1+1:end) > this_right,1,'first');
          if isempty(ii), ii = length(t_eye)-idx1; end
          idx2 = idx1 + ii;
          
          if idx2 > idx1

              % Only take window around initiation of pass
%               idx_s = max(1,idx1+window(1));
%               idx_e = min(idx2, idx1+window(2));
%               idx1 = idx_s; idx2 = idx_e;
              
              idx_passing(idx1:idx2) = true;
              passing_intervals(end+1,:) = [idx1 idx2];
              passing_times(end+1,:) = [this_left this_right];
              passing_cycles(idx1:idx2) = c;

              % Assign difficulty to this interval
              passing_diff(idx1:idx2) = get_difficulty(seq_diff.difficulty, preprocess.sim2track.matrix, ...
                                                       mean([this_left this_right]), ...
                                                       params.lane_dist);

              % Assign outcome to this interval
              outcome = get_outcome(preprocess.sim2track, ...
                                                       [this_left this_right], ...
                                                       params.lane_dist);

              passing_outcomes(idx1:idx2) = repmat(outcome,idx2-idx1+1,1);

              results.overtake_outcomes(end+1) = outcome;
          
          end
                                               
      end
      % Reset
      this_left = -1;
      this_right = -1;
   end
end

% Compute mean/variance for all data
% if isfield(preprocess, 'luminance') && ~preprocess.luminance.deficient
%     pd = preprocess.luminance.diam_left{1};
% else
%     pd = data.eye.diam_left;
% end
pdz = zscore(pd);

sr = preprocess.saccades.saccade_rate;

results.diff_levels = diff_levels;

results.baseline.pupil = pd(idx_baseline);
results.nobaseline.pupil = pd(~idx_baseline);
results.baseline.saccade_rate = sr(idx_baseline);
results.nobaseline.saccade_rate = sr(~idx_baseline);

results.passing.pupil = pd(idx_passing);
results.passing.saccade_rate = sr(idx_passing);

N_cycles = length(preprocess.sim2track.cycle_times)+1;
results.cycles.baseline.pupil = cell(N_cycles,1);
for i = 1 : N_cycles
    results.cycles.baseline.pupil(i) = {pd(idx_baseline & baseline_cycles==i)};
    results.cycles.nobaseline.pupil(i) = {pd(~idx_baseline & baseline_cycles==i)};
    results.cycles.passing.pupil(i) = {pd(idx_passing & passing_cycles==i)};
    
    results.zscore.cycles.baseline.pupil(i) = {pdz(idx_baseline & baseline_cycles==i)};
    results.zscore.cycles.nobaseline.pupil(i) = {pdz(~idx_baseline & baseline_cycles==i)};
    results.zscore.cycles.passing.pupil(i) = {pdz(idx_passing & passing_cycles==i)};
end

results.zscore.baseline.pupil = pdz(idx_baseline);
results.zscore.nobaseline.pupil = pdz(~idx_baseline);

results.zscore.passing.pupil = pdz(idx_passing);

results.passing_diff.pupil = cell(length(diff_levels),1);
results.zscore.passing_diff.pupil = cell(length(diff_levels),1);

results.passing_diff.cycles.passing.pupil = cell(length(diff_levels), N_cycles);
results.zscore.passing_diff.cycles.passing.pupil = cell(length(diff_levels), N_cycles);

for i = 1 : length(diff_levels)
    results.passing_diff.pupil(i) = {pd(passing_diff==diff_levels(i))};
    results.zscore.passing_diff.pupil(i) = {pdz(passing_diff==diff_levels(i))};
    
    for j = 1 : N_cycles
        results.passing_diff.cycles.passing.pupil(i,j) = {pd(idx_passing & passing_cycles==j & passing_diff==diff_levels(i))};
        results.zscore.passing_diff.cycles.passing.pupil(i,j) = {pdz(idx_passing & passing_cycles==j & passing_diff==diff_levels(i))};
    end
    
end

results.idx_baseline = idx_baseline;
results.idx_passing = idx_passing;
results.passing_diffs = passing_diff;
results.passing_outcomes = passing_outcomes;

% Outcomes
results.passing_outcome.positive.pupil = pd(passing_outcomes>0);
results.zscore.passing_outcome.positive.pupil = pdz(passing_outcomes>0);
results.passing_outcome.negative.pupil = pd(passing_outcomes<0);
results.zscore.passing_outcome.negative.pupil = pdz(passing_outcomes<0);

% Add each baseline interval
results.intervals.baseline.idx = baseline_intervals;
results.intervals.baseline.times = zeros(size(baseline_intervals,1),1);
results.intervals.baseline.pupil = cell(size(baseline_intervals,1),1);
results.intervals.baseline.saccade_rate = cell(size(baseline_intervals,1),1);
results.zscore.intervals.baseline.pupil = cell(size(baseline_intervals,1),1);
for i = 1 : size(baseline_intervals,1)
    results.intervals.baseline.times(i) = mean(baseline(i,:));
    results.intervals.baseline.pupil(i) = ...
        {pd(baseline_intervals(i,1):baseline_intervals(i,2))};    
    results.zscore.intervals.baseline.pupil(i) = ...
        {pdz(baseline_intervals(i,1):baseline_intervals(i,2))}; 
    results.intervals.baseline.saccade_rate(i) = ...
        {sr(baseline_intervals(i,1):baseline_intervals(i,2))};
end

% Add each passing interval
results.intervals.passing.idx = passing_intervals;
results.intervals.passing.times = zeros(size(passing_intervals,1),1);
results.intervals.passing.pupil = cell(size(passing_intervals,1),1);
results.intervals.passing.saccade_rate = cell(size(passing_intervals,1),1);
results.zscore.intervals.passing.pupil = cell(size(passing_intervals,1),1);
for i = 1 : size(passing_intervals,1)
    results.intervals.passing.times(i) = mean(passing_times(i,:));
    results.intervals.passing.pupil(i) = ...
        {pd(passing_intervals(i,1):passing_intervals(i,2))};  
    results.zscore.intervals.passing.pupil(i) = ...
        {pdz(passing_intervals(i,1):passing_intervals(i,2))};  
    results.intervals.passing.saccade_rate(i) = ...
        {sr(passing_intervals(i,1):passing_intervals(i,2))};  
end

% Compute mean/variance for each cycle/round 



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
