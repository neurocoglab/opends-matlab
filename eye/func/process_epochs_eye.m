function [ results, summary ] = process_epochs_eye( params, data, summary )
% Summarize data over epochs
%
% Arguments:
%
% data:     output from preprocessing_eye
%
% Epochs:   baseline, non-baseline, passing events
% Data:     pupil diameter, saccade rate, blink rate
%

hdr_epochs = {'SubjectID','EpochID','EpochType','TimeStart','TimeEnd','Cycle','MeanPupil','Outcome'};

if params.sim.epochs.difficulty.apply
    hdr_epochs = [hdr_epochs;{'Difficulty'}];
end

if isempty(summary)
    summary.subjects = [];
    summary.baseline.pupil = [];
    summary.passing.pupil  = [];
    summary.zscore.baseline.pupil = [];
    summary.zscore.passing.pupil  = [];
    summary.passing_diff.pupil  = [];
    summary.zscore.passing_diff.pupil  = [];
    summary.baseline.saccade_rate = [];
    summary.passing.saccade_rate  = [];
    summary.cycles.baseline.pupil = [];
    summary.zscore.cycles.baseline.pupil = [];
    summary.cycles.passing.pupil = [];
    summary.zscore.cycles.passing.pupil = [];
    summary.passing_outcome.positive.pupil = [];
    summary.passing_outcome.negative.pupil = [];
    summary.zscore.passing_outcome.positive.pupil = [];
    summary.zscore.passing_outcome.negative.pupil = [];
    summary.passing_diffs = {};
    summary.passing_outcomes = {};
    
    summary.idx_baseline = [];
    summary.idx_passing = [];
    
    summary.baseline.subjects.pupil = [];
    summary.passing.subjects.pupil = [];
    summary.zscore.baseline.subjects.pupil = [];
    summary.zscore.passing.subjects.pupil = [];

    summary.cycles.baseline.subjects.pupil = [];
    summary.cycles.passing.subjects.pupil = [];
    summary.zscore.cycles.baseline.subjects.pupil = [];
    summary.zscore.cycles.passing.subjects.pupil = [];
    
    summary.passing_outcome.positive.subjects.pupil = [];
    summary.passing_outcome.negative.subjects.pupil = [];
    summary.zscore.passing_outcome.positive.subjects.pupil = [];
    summary.zscore.passing_outcome.negative.subjects.pupil = [];

    summary.baseline.subjects.saccade_rate = [];
    summary.passing.subjects.saccade_rate = [];

    summary.epochs_table = [];
    
end

results = [];
results.subject = data.subject;

epochs_table = cell2table(cell(0,length(hdr_epochs)), 'VariableNames', hdr_epochs);
current_epoch = 1;

% Load sequence difficulty ratings
if params.sim.epochs.difficulty.apply
    seq_diff = readtable(sprintf('%s/%s/%s', params.io.input_dir, ...
                                             params.io.metadata_dir, ...
                                             params.sim.epochs.difficulty.sequence_file));
    diff_levels = params.sim.epochs.difficulty.levels; % [1 2]; % unique(D);
end

% Baseline
baseline = data.sim.sim2track.baseline;
 
if isfield(data.sim, 'luminance') && ~data.sim.luminance.deficient
    idx_offset = params.eye.luminance.use_offset;
    if idx_offset < 1
       [~,idx_offset] = max(data.eye.luminance.r2); 
    end
    t_eye = data.eye.luminance.ts{idx_offset};
    pd = data.eye.luminance.diam{idx_offset};
else
    t_eye = data.eye.t;
    if isfield(data.eye, 'blinks')
        pd = data.eye.blinks.diam;
    else
        pd = data.eye.diam;
    end
end

pd = pd(:);
idx_nan = isnan(pd);
pdz = pd;
pdz(idx_nan) = mean(pdz(~idx_nan));
pdz = zscore(pdz);
pdz(idx_nan) = nan;

sr = data.eye.saccades.saccade_rate;

N = length(t_eye);
idx_baseline = false(N,1);
baseline_intervals = zeros(0,2);
baseline_cycles = zeros(N,1);

for i = 1 : height(baseline)
   
    % Map baseline times to ts indexes
    ti = baseline.Start(i);
    c = find(data.sim.sim2track.cycle_times > ti, 1);
    if isempty(c); c = length(data.sim.sim2track.cycle_times)+1; end
    idx1 = find(t_eye < ti,1,'last');
    idx2 = idx1 + find(t_eye(idx1+1:end) > baseline.End(i),1,'first');
    
    if idx2 > idx1
        idx_baseline(idx1:idx2) = true;
        baseline_cycles(idx1:idx2) = c;
        try
        baseline_intervals(end+1,:) = [idx1 idx2];
        catch
           a=0; 
        end
        % Append to table
        cycle = get_cycle(data.sim.sim2track, data.eye.t(idx1));
        row_data = {data.subject, current_epoch, 'Baseline', data.eye.t(idx1), data.eye.t(idx2), cycle, mean(pdz(idx1:idx2)), ''};
        if params.sim.epochs.difficulty.apply
            row_data = [row_data {''}];
        end
        row = cell2table(row_data, 'VariableNames', hdr_epochs);
        epochs_table = [epochs_table;row];
        current_epoch = current_epoch + 1;
    end
    
end

results.eye.epochs.idx_baseline = idx_baseline;

% Passing
passing_intervals = zeros(0,2);
passing_times = zeros(0,2);
idx_passing = false(N,1);
passing_diff = zeros(N,1);
passing_outcomes = zeros(N,1);
passing_cycles = zeros(N,1);

max_interval = 30000; % Half a minute
        
left = [data.sim.sim2track.left_change_times, ...
        true(length(data.sim.sim2track.left_change_times),1)];
right = [data.sim.sim2track.right_change_times, ...
         false(length(data.sim.sim2track.right_change_times),1)];

left_right = [left;right];
[~,idx] = sort(left_right(:,1));
left_right = left_right(idx,:);

results.eye.epochs.overtake_intervals = zeros(0,2);
results.eye.epochs.overtake_outcomes = [];

this_left = -1;
this_right = -1;
window = round(params.eye.epochs.overtake_window * data.eye.Fs);
for j = 1 : length(left_right)
   if left_right(j,2)
      % Is change to left 
      this_left = left_right(j,1);
   else
      % Is change to right
      this_right = left_right(j,1);
      if this_left > 0 && this_right-this_left < max_interval
          % Valid passing segment, add
          results.eye.epochs.overtake_intervals(end+1,:) = [this_left this_right];
          
          c = find(data.sim.sim2track.cycle_times > this_left, 1);
          if isempty(c); c = length(data.sim.sim2track.cycle_times)+1; end
          
          % Map lane change times to ts indexes
          idx1 = find(t_eye < this_left,1,'last');
          ii = find(t_eye(idx1+1:end) > this_right,1,'first');
          if isempty(ii), ii = length(t_eye)-idx1; end
          idx2 = idx1 + ii;
          
          if idx2 > idx1

              % Only take window around initiation of pass
              idx_passing(idx1:idx2) = true;
              passing_intervals(end+1,:) = [idx1 idx2];
              passing_times(end+1,:) = [this_left this_right];
              passing_cycles(idx1:idx2) = c;

              % Assign difficulty to this interval
              if params.sim.epochs.difficulty.apply
                  passing_diff(idx1:idx2) = get_difficulty(seq_diff, data.sim.sim2track.matrix, ...
                                                           mean([this_left this_right]), ...
                                                           params.sim.lane_dist);
              end
              
              % Assign outcome to this interval
              outcome = get_outcome(data.sim.sim2track, [this_left this_right]);
              cycle = get_cycle(data.sim.sim2track, this_left);

              passing_outcomes(idx1:idx2) = repmat(outcome,idx2-idx1+1,1);

              results.eye.epochs.overtake_outcomes(end+1) = outcome;
          
              % Append to table
              row_data = {data.subject, current_epoch, 'Passing', data.eye.t(idx1), data.eye.t(idx2), cycle, mean(pdz(idx1:idx2)), outcome};
              if params.sim.epochs.difficulty.apply
                  row_data = [row_data {passing_diff(idx1)}];
              end
              row = cell2table(row_data, 'VariableNames', hdr_epochs);
              epochs_table = [epochs_table;row];
              current_epoch = current_epoch + 1;

          end
                                               
      end
      % Reset
      this_left = -1;
      this_right = -1;
   end
end

% Compute mean/variance for all data
% if isfield(preprocess, 'luminance') && ~data.sim.luminance.deficient
%     pd = data.sim.luminance.diam_left{1};
% else
%     pd = data.eye.diam_left;
% end


results.eye.epochs.baseline.pupil = pd(idx_baseline);
results.eye.epochs.nobaseline.pupil = pd(~idx_baseline);
results.eye.epochs.baseline.saccade_rate = sr(idx_baseline);
results.eye.epochs.nobaseline.saccade_rate = sr(~idx_baseline);

results.eye.epochs.passing.pupil = pd(idx_passing);
results.eye.epochs.passing.saccade_rate = sr(idx_passing);

N_cycles = length(data.sim.sim2track.cycle_times)+1;
results.eye.epochs.cycles.baseline.pupil = cell(N_cycles,1);
for i = 1 : N_cycles
    results.eye.epochs.cycles.baseline.pupil(i) = {pd(idx_baseline & baseline_cycles==i)};
    results.eye.epochs.cycles.nobaseline.pupil(i) = {pd(~idx_baseline & baseline_cycles==i)};
    results.eye.epochs.cycles.passing.pupil(i) = {pd(idx_passing & passing_cycles==i)};
    
    results.eye.epochs.zscore.cycles.baseline.pupil(i) = {pdz(idx_baseline & baseline_cycles==i)};
    results.eye.epochs.zscore.cycles.nobaseline.pupil(i) = {pdz(~idx_baseline & baseline_cycles==i)};
    results.eye.epochs.zscore.cycles.passing.pupil(i) = {pdz(idx_passing & passing_cycles==i)};
end

results.eye.epochs.zscore.baseline.pupil = pdz(idx_baseline);
results.eye.epochs.zscore.nobaseline.pupil = pdz(~idx_baseline);

results.eye.epochs.zscore.passing.pupil = pdz(idx_passing);

if params.sim.epochs.difficulty.apply
    results.eye.epochs.diff_levels = diff_levels;
    
    results.eye.epochs.passing_diff.pupil = cell(length(diff_levels),1);
    results.eye.epochs.zscore.passing_diff.pupil = cell(length(diff_levels),1);

    results.eye.epochs.passing_diff.cycles.passing.pupil = cell(length(diff_levels), N_cycles);
    results.eye.epochs.zscore.passing_diff.cycles.passing.pupil = cell(length(diff_levels), N_cycles);

    for i = 1 : length(diff_levels)
        results.eye.epochs.passing_diff.pupil(i) = {pd(passing_diff==diff_levels(i))};
        results.eye.epochs.zscore.passing_diff.pupil(i) = {pdz(passing_diff==diff_levels(i))};

        for j = 1 : N_cycles
            results.eye.epochs.passing_diff.cycles.passing.pupil(i,j) = {pd(idx_passing & passing_cycles==j & passing_diff==diff_levels(i))};
            results.eye.epochs.zscore.passing_diff.cycles.passing.pupil(i,j) = {pdz(idx_passing & passing_cycles==j & passing_diff==diff_levels(i))};
        end

    end
    results.eye.epochs.passing_diffs = passing_diff;
end

results.eye.epochs.idx_baseline = idx_baseline;
results.eye.epochs.idx_passing = idx_passing;
results.eye.epochs.passing_outcomes = passing_outcomes;

% Outcomes
results.eye.epochs.passing_outcome.positive.pupil = pd(passing_outcomes>0);
results.eye.epochs.zscore.passing_outcome.positive.pupil = pdz(passing_outcomes>0);
results.eye.epochs.passing_outcome.negative.pupil = pd(passing_outcomes<0);
results.eye.epochs.zscore.passing_outcome.negative.pupil = pdz(passing_outcomes<0);

% Add each baseline interval
results.eye.epochs.intervals.baseline.idx = baseline_intervals;
results.eye.epochs.intervals.baseline.times = zeros(size(baseline_intervals,1),1);
results.eye.epochs.intervals.baseline.pupil = cell(size(baseline_intervals,1),1);
results.eye.epochs.intervals.baseline.saccade_rate = cell(size(baseline_intervals,1),1);
results.eye.epochs.zscore.intervals.baseline.pupil = cell(size(baseline_intervals,1),1);
for i = 1 : size(baseline_intervals,1)
    results.eye.epochs.intervals.baseline.times(i) = mean([baseline.Start(i),baseline.End(i)]);
    results.eye.epochs.intervals.baseline.pupil(i) = ...
        {pd(baseline_intervals(i,1):baseline_intervals(i,2))};    
    results.eye.epochs.zscore.intervals.baseline.pupil(i) = ...
        {pdz(baseline_intervals(i,1):baseline_intervals(i,2))}; 
    results.eye.epochs.intervals.baseline.saccade_rate(i) = ...
        {sr(baseline_intervals(i,1):baseline_intervals(i,2))};
end

% Add each passing interval
results.eye.epochs.intervals.passing.idx = passing_intervals;
results.eye.epochs.intervals.passing.times = zeros(size(passing_intervals,1),1);
results.eye.epochs.intervals.passing.pupil = cell(size(passing_intervals,1),1);
results.eye.epochs.intervals.passing.saccade_rate = cell(size(passing_intervals,1),1);
results.eye.epochs.zscore.intervals.passing.pupil = cell(size(passing_intervals,1),1);
for i = 1 : size(passing_intervals,1)
    results.eye.epochs.intervals.passing.times(i) = mean(passing_times(i,:));
    results.eye.epochs.intervals.passing.pupil(i) = ...
        {pd(passing_intervals(i,1):passing_intervals(i,2))};  
    results.eye.epochs.zscore.intervals.passing.pupil(i) = ...
        {pdz(passing_intervals(i,1):passing_intervals(i,2))};  
    results.eye.epochs.intervals.passing.saccade_rate(i) = ...
        {sr(passing_intervals(i,1):passing_intervals(i,2))};  
end

% Compute mean/variance for each cycle/round 

 % Aggregate summary stats
summary.subjects = [summary.subjects {data.subject}];

summary.epochs_table = [summary.epochs_table;epochs_table];

summary.baseline.subjects.pupil = [summary.baseline.subjects.pupil {results.eye.epochs.baseline.pupil}];
summary.passing.subjects.pupil = [summary.passing.subjects.pupil {results.eye.epochs.passing.pupil}];
summary.zscore.baseline.subjects.pupil = [summary.zscore.baseline.subjects.pupil {results.eye.epochs.zscore.baseline.pupil}];
summary.zscore.passing.subjects.pupil = [summary.zscore.passing.subjects.pupil {results.eye.epochs.zscore.passing.pupil}];

summary.cycles.baseline.subjects.pupil = [summary.cycles.baseline.subjects.pupil {results.eye.epochs.cycles.baseline.pupil}];
summary.cycles.passing.subjects.pupil = [summary.cycles.passing.subjects.pupil {results.eye.epochs.cycles.passing.pupil}];
summary.zscore.cycles.baseline.subjects.pupil = [summary.zscore.cycles.baseline.subjects.pupil {results.eye.epochs.zscore.cycles.baseline.pupil}];
summary.zscore.cycles.passing.subjects.pupil = [summary.zscore.cycles.passing.subjects.pupil {results.eye.epochs.zscore.cycles.passing.pupil}];

summary.baseline.subjects.saccade_rate = [summary.baseline.subjects.saccade_rate {results.eye.epochs.baseline.saccade_rate}];
summary.passing.subjects.saccade_rate = [summary.passing.subjects.saccade_rate {results.eye.epochs.passing.saccade_rate}];

N_cycles = params.sim.rounds.max_cycles; % length(data.sim.sim2track.cycle_times)+1;
N_cycles_subj = length(data.sim.sim2track.cycle_times)+1;

if isempty(summary.cycles.baseline.pupil)
    summary.cycles.baseline.pupil = cell(N_cycles,1);
    summary.zscore.cycles.baseline.pupil = cell(N_cycles,1);
    summary.cycles.passing.pupil = cell(N_cycles,1);
    summary.zscore.cycles.passing.pupil = cell(N_cycles,1);

end

for j = 1 : N_cycles
    if j <= N_cycles_subj
        summary.cycles.baseline.pupil(j) = {[summary.cycles.baseline.pupil{j};results.eye.epochs.cycles.baseline.pupil{j}]};
        summary.zscore.cycles.baseline.pupil(j) = {[summary.zscore.cycles.baseline.pupil{j};results.eye.epochs.zscore.cycles.baseline.pupil{j}]};
        summary.cycles.passing.pupil(j) = {[summary.cycles.passing.pupil{j};results.eye.epochs.cycles.passing.pupil{j}]};
        summary.zscore.cycles.passing.pupil(j) = {[summary.zscore.cycles.passing.pupil{j};results.eye.epochs.zscore.cycles.passing.pupil{j}]};
    else
        summary.cycles.baseline.pupil(j) = {[summary.cycles.baseline.pupil{j};{}]};
        summary.zscore.cycles.baseline.pupil(j) = {[summary.zscore.cycles.baseline.pupil{j};{}]};
        summary.cycles.passing.pupil(j) = {[summary.cycles.passing.pupil{j};{}]};
        summary.zscore.cycles.passing.pupil(j);
    end
end

if params.sim.epochs.difficulty.apply
    summary.diff_levels = results.eye.epochs.diff_levels;
    
    if isempty(summary.passing_diff.pupil)
        summary.passing_diff.pupil = results.eye.epochs.passing_diff.pupil;
        summary.zscore.passing_diff.pupil = results.eye.epochs.zscore.passing_diff.pupil;
    else
        for j = 1 : length(results.eye.epochs.diff_levels)
            summary.passing_diff.pupil(j) = {[summary.passing_diff.pupil{j};results.eye.epochs.passing_diff.pupil{j}]};
            summary.zscore.passing_diff.pupil(j) = {[summary.zscore.passing_diff.pupil{j};results.eye.epochs.zscore.passing_diff.pupil{j}]};
        end
    end

    summary.passing_diffs = [summary.passing_diffs {results.eye.epochs.passing_diffs}];

end

summary.idx_baseline = [summary.idx_baseline {find(results.eye.epochs.idx_baseline)}];
summary.idx_passing = [summary.idx_passing {find(results.eye.epochs.idx_passing)}];

summary.passing_outcomes = [summary.passing_outcomes {results.eye.epochs.passing_outcomes}];
summary.passing_outcome.positive.subjects.pupil = [summary.passing_outcome.positive.subjects.pupil {results.eye.epochs.passing_outcome.positive.pupil}];
summary.passing_outcome.negative.subjects.pupil = [summary.passing_outcome.negative.subjects.pupil {results.eye.epochs.passing_outcome.negative.pupil}];
summary.zscore.passing_outcome.positive.subjects.pupil = [summary.zscore.passing_outcome.positive.subjects.pupil {results.eye.epochs.zscore.passing_outcome.positive.pupil}];
summary.zscore.passing_outcome.negative.subjects.pupil = [summary.zscore.passing_outcome.negative.subjects.pupil {results.eye.epochs.zscore.passing_outcome.negative.pupil}];

summary.baseline.pupil = [summary.baseline.pupil;results.eye.epochs.baseline.pupil];
summary.passing.pupil = [summary.passing.pupil;results.eye.epochs.passing.pupil];

summary.baseline.saccade_rate = [summary.baseline.saccade_rate;results.eye.epochs.baseline.saccade_rate];
summary.passing.saccade_rate = [summary.passing.saccade_rate;results.eye.epochs.passing.saccade_rate];


end

