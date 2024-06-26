function [ results ] = process_hilbert_eeg ( params, data, results )
%%%%%%%%%%%%%%
% Compute Hilbert transforms (amplitude envelopes) for specified frequency bands
% from the given EEG data
%

if ~exist('results', 'var') 
    results = [];
end

results.subject = data.subject;

outdir = sprintf( '%s/%s', params.io.output_dir, data.subject );
figdir = sprintf( '%s/figures', outdir );

cfg_hb = data.eeg.cfg;
cfg_hb.continuous = 1;
cfg_hb.demean = 'no';

% Load Hilbert bands and filter orders from CSV file
fn = params.eeg.hilbert.bands_file;
T = readtable( fn );

N_bands = height(T);
N_channels = length(data.eeg.all_channels);
channel_labels = sort(data.eeg.all_channels);
N_sim = length(data.sim.sim2track.cycle_times);
N_cycles = min(params.sim.rounds.max_cycles, N_sim+1);
t_eeg = data.eeg.ft.time{1} * 1000; % Time is milliseconds
t_max = t_eeg(end);
if N_cycles < params.sim.rounds.max_cycles
    warning('Subject has fewer than max cycles (%d<%d)', N_cycles, params.sim.rounds.max_cycles)
elseif N_sim >= params.sim.rounds.max_cycles
    warning('Subject has more than max cycles (%d>%d); will truncate to max', ...
                N_sim+1, params.sim.rounds.max_cycles)
    t_max = data.sim.sim2track.cycle_times(N_cycles);
end
idx_keep = t_eeg <= t_max;
t_eeg = t_eeg(idx_keep);

results.eeg.hilbert.bands = T;
results.eeg.hilbert.epochs.baselines = cell(N_bands,1);
results.eeg.hilbert.epochs.overtakes = cell(N_bands,1);
results.eeg.hilbert.N_cycles = N_cycles;
results.eeg.hilbert.epochs.stats.baseline.mean = nan(N_bands,N_channels);
results.eeg.hilbert.epochs.stats.overtake.mean = nan(N_bands,N_channels);
results.eeg.hilbert.epochs.stats.cycles.mean = nan(N_bands,N_channels,N_cycles);
results.eeg.hilbert.epochs.stats.baseline.cycles.mean = nan(N_bands,N_channels,N_cycles);
results.eeg.hilbert.epochs.stats.overtake.cycles.mean = nan(N_bands,N_channels,N_cycles);
if params.sim.epochs.difficulty.apply
    results.eeg.hilbert.epochs.stats.overtake.difficulty.mean = ...
        nan(N_bands,N_channels,length(params.sim.epochs.difficulty.levels));
end
if params.sim.epochs.outcomes.apply
    results.eeg.hilbert.epochs.stats.overtake.outcomes.mean = ...
        nan(N_bands,N_channels,length(params.sim.epochs.outcomes.levels));
end
results.eeg.hilbert.envelopes = [];
results.eeg.hilbert.envelopes.raw = cell(N_bands,1);
results.eeg.hilbert.envelopes.zscore = cell(N_bands,1);
results.eeg.hilbert.envelopes.time = [];

idata = data.eeg.ft;
Fs = idata.fsample;

% Load sequence difficulty ratings
if params.sim.epochs.difficulty.apply
    seq_diff = readtable(sprintf('%s/%s/%s', params.io.input_dir, ...
                                             params.io.metadata_dir, ...
                                             params.sim.epochs.difficulty.sequence_file));
    diff_levels = params.sim.epochs.difficulty.levels; % [1 2]; % unique(D);
end

for bb = 1 : N_bands

    freqband = [T.From(bb) T.To(bb)];
    fprintf('\tComputing envelopes for %s: [%1.1f to %1.1f Hz]', ...
                T.Band{bb}, freqband(1), freqband(2));

    % make filter
    cfg_hb = [];
    cfg_hb.filtord = T.Filtord(bb);
    [b,a] = butter(cfg_hb.filtord,2*freqband/idata.fsample, 'bandpass');

    if params.eeg.hilbert.plots.save
        fmin = max(0.1,freqband(1)-10);
        fmax = max(20, freqband(2)+10);
        freqz(b,a,fmin:(fmax-fmin)/1000:fmax,idata.fsample);
        saveas(gcf,sprintf('%s/eeg_hilbert_filter_%s.png', figdir, T.Band{bb}));
        close(gcf);
    end

    % Apply filter to data
    eeg_channels = data.eeg.eeg_channels;
    Fs2 = Fs;

    if params.eeg.hilbert.save_filtered
        datfilt = cell(N_channels,1);
    end
    envelopes = cell(N_channels,1);
    envelopes_zscore = cell(N_channels,1);
    %times = cell(N_channels,1);
    eeg_channels_sorted = cell(N_channels,1);
    
    for cc = 1 : N_channels
        channel = channel_labels{cc};
        idx = find(strcmp(idata.label, channel));
        if isempty(idx)
            %eeg_channels_sorted(cc)={};
            continue;
        end
        eeg_channels_sorted(cc)={channel};
        X = idata.trial{1}(idx,idx_keep);
        idx_nan = ~isnan(X);
        if isempty(results.eeg.hilbert.envelopes.time)
            t_eeg = t_eeg(idx_nan);
        end
        Y = filtfilt(b,a,X(idx_nan));
        
        datfilt_cc = Y;
        H = hilbert(Y);
        
         % Downsample?
        if params.eeg.hilbert.downsample > 0 && params.eeg.hilbert.downsample < Fs
            step = round(idata.fsample/params.eeg.hilbert.downsample);
            ts_env = timeseries(H, t_eeg);
            idx_ds = 1:step:length(ts_env.Time);
            ts_env = resample(ts_env, ts_env.Time(idx_ds));
            H = ts_env.Data;
            if params.eeg.hilbert.save_filtered
                ts_filt = timeseries(datfilt_cc, t_eeg);
                ts_filt = resample(ts_filt, ts_filt.Time(idx_ds));
                datfilt_cc = ts_filt.Data;
            end
            Fs2 = 1/(ts_env.Time(2)- ts_env.Time(1)); % Time is in seconds
            t_eeg_ds = single(ts_env.Time);
        else
            t_eeg_ds = t_eeg;
        end
        
        if isempty(results.eeg.hilbert.envelopes.time)
            results.eeg.hilbert.envelopes.time = t_eeg_ds; % Assuming that nans are not channel-specific
        end

        envelopes(cc) = {squeeze(single(H))};
        envelopes_zscore(cc) = {zscore(abs(envelopes{cc}))};
        %times(cc) = {time_cc};
        if params.eeg.hilbert.save_filtered
            datfilt(cc) = {single(datfilt_cc)};
        end
    end
    
    if params.general.debug
        fprintf(' [downsampled from %1.2f to %1.2f Hz]', Fs, Fs2);
    end
            
    t_eeg_ds = results.eeg.hilbert.envelopes.time;
    results.eeg.hilbert.cfg(bb) = cfg_hb;
    results.eeg.hilbert.Fs = Fs2;
    if params.eeg.hilbert.save_filtered
        results.eeg.hilbert.filtered(bb) = {datfilt};
    end
    results.eeg.hilbert.channels = eeg_channels_sorted;
    results.eeg.hilbert.envelopes.raw(bb) = {envelopes};
    results.eeg.hilbert.envelopes.zscore(bb) = {envelopes_zscore};
    %results.eeg.hilbert.time(bb) = {times};

    % Baseline vs. overtakes 
    baseline_idx = cell(N_channels,1);
    overtake_idx = cell(N_channels,1);
    for cc = 1 : N_channels
        channel = eeg_channels_sorted{cc};
        if isempty(channel)
            % Will be nan's for this missing channel
            continue;
        end

        %t_eeg = times{cc} * 1000; % Seconds to milliseconds
        
        % Map baseline times to EEG time series indices
        baseline_info = get_baseline_info(data, t_eeg_ds, N_cycles, cc==1);
        baseline_idx{cc} = baseline_info.indices;
        results.eeg.hilbert.epochs.stats.baseline.mean(bb,cc) = mean(envelopes_zscore{cc}(baseline_info.indices));

        % Map overtake times to EEG time series indexes
        overtake_info = get_overtake_info(data, params, t_eeg_ds);
        overtake_idx{cc} = overtake_info.indices;
        results.eeg.hilbert.epochs.stats.overtake.mean(bb,cc) = mean(envelopes_zscore{cc}(overtake_info.indices));

        % Overtake - difficulty
        if params.sim.epochs.difficulty.apply
            results.eeg.hilbert.epochs.overtake.difficulty = overtake_info.difficulty;
            for k = 1 : length(params.sim.epochs.difficulty.levels)
                level = params.sim.epochs.difficulty.levels(k);
                results.eeg.hilbert.epochs.stats.overtake.difficulty.mean(bb,cc,k) = ...
                    mean(envelopes_zscore{cc}(overtake_info.difficulty==level));
            end
        end

        % Overtake - outcomes
        if params.sim.epochs.outcomes.apply
            results.eeg.hilbert.epochs.overtake.outcome = overtake_info.outcomes;
            for k = 1 : length(params.sim.epochs.outcomes.levels)
                level = params.sim.epochs.outcomes.levels(k);
                results.eeg.hilbert.epochs.stats.overtake.outcomes.mean(bb,cc,k) = ...
                    mean(envelopes_zscore{cc}(overtake_info.outcomes==level));
            end
        end

        % Stats per cycle
        cycles = get_cycles(data, t_eeg_ds, N_cycles);
        for c = 1 : N_cycles
            results.eeg.hilbert.epochs.stats.cycles.mean(bb,cc,c) = mean(envelopes_zscore{cc}(cycles==c));
            results.eeg.hilbert.epochs.stats.baseline.cycles.mean(bb,cc,c) = mean(envelopes_zscore{cc}(baseline_info.cycles==c));
            results.eeg.hilbert.epochs.stats.overtake.cycles.mean(bb,cc,c) = mean(envelopes_zscore{cc}(overtake_info.cycles==c));
        end

    end
    
    results.eeg.hilbert.epochs.baselines{bb} = baseline_idx;
    results.eeg.hilbert.epochs.overtakes{bb} = overtake_idx;
    
    fprintf('\n');

end

end

% Extract baseline information in terms of t_eeg
function baseline_info = get_baseline_info(data, t_eeg, N_cycles, warn)

    if nargin < 4
        warn = false;
    end
    baseline = data.sim.sim2track.baseline;
    N = length(t_eeg);
    baseline_info = [];
    baseline_info.indices = false(N,1);
    baseline_info.intervals = zeros(0,2);
    baseline_info.cycles = zeros(N,1);
    t_last = 0;
    for i = 1 : N_cycles %size(baseline,1)
        baseline_i = baseline(baseline.Cycle==i,:);
        if height(baseline_i)==0
            if warn
                warning('Missing baseline for cycle %d; check log?', i);
            end
            continue;
        end
        % Take first baseline period regardless of repeats
        ti = baseline_i.Start(1);
        tj = baseline_i.End(1);
        if ti < t_last
            % Weird, skip
            if warn
                warning('Skipping overlapping baseline interval (cycle %d); check log?', i);
            end
            continue;
        end
        c = find(data.sim.sim2track.cycle_times > ti, 1);
        if isempty(c); c = N_cycles; end
        idx1 = find(t_eeg < ti, 1, 'last');
        idx2 = idx1 + find(t_eeg(idx1+1:end) > tj, 1,'first');
        
        if idx2 > idx1
            baseline_info.indices(idx1:idx2) = true;
            baseline_info.cycles(idx1:idx2) = c;
            baseline_info.intervals(end+1,:) = [idx1 idx2];
        end

        t_last = tj;
    
    end

end

function cycles = get_cycles(data, t_eeg, N_cycles)

    N_sim = length(data.sim.sim2track.cycle_times);
    N_t = length(t_eeg);
    cycles = nan(N_t,1);
    idx1 = 1;
    for i = 1 : N_cycles
        if i <= N_sim
            cti = data.sim.sim2track.cycle_times(i);
            idx2 = find(t_eeg > cti, 1);
            if isempty(idx2), idx2 = length(t_eeg); end
            cycles(idx1:idx2-1) = i;
            idx1 = idx2; 
        end
    end
    if idx2 < N_t && N_cycles > N_sim
        cycles(idx2:end) = N_cycles;
    end

end

% Extract overtake information in terms of t_eeg
function overtake_info = get_overtake_info(data, params, t_eeg)

    N = length(t_eeg);

    overtake_info = [];
    overtake_info.intervals = zeros(0,2);
    overtake_info.passing_times = zeros(0,2);
    overtake_info.indices = false(N,1);
    overtake_info.difficulty = zeros(N,1);
    overtake_info.outcomes = zeros(N,1);
    overtake_info.cycles = zeros(N,1);
    
    max_interval = 30000; % Half a minute
            
    left = [data.sim.sim2track.left_change_times, ...
            true(length(data.sim.sim2track.left_change_times),1)];
    right = [data.sim.sim2track.right_change_times, ...
             false(length(data.sim.sim2track.right_change_times),1)];
    
    left_right = [left;right];
    [~,idx] = sort(left_right(:,1));
    left_right = left_right(idx,:);
    
    this_left = -1;
    for j = 1 : length(left_right)
       if left_right(j,2)
          % Is change to left 
          this_left = left_right(j,1);
       else
          % Is change to right
          this_right = left_right(j,1);
          if this_left > 0 && this_right-this_left < max_interval
              % Valid passing segment, add
              overtake_info.intervals(end+1,:) = [this_left this_right];
              
              c = find(data.sim.sim2track.cycle_times > this_left, 1);
              if isempty(c); c = length(data.sim.sim2track.cycle_times)+1; end
              
              % Map lane change times to ts indexes
              idx1 = find(t_eeg < this_left,1,'last');
              ii = find(t_eeg(idx1+1:end) > this_right,1,'first');
              if isempty(ii), ii = length(t_eeg)-idx1; end
              idx2 = idx1 + ii;
              
              if idx2 > idx1
    
                  % Only take window around initiation of pass
                  overtake_info.indices(idx1:idx2) = true;
                  overtake_info.intervals(end+1,:) = [idx1 idx2];
                  overtake_info.passing_times(end+1,:) = [this_left this_right];
                  overtake_info.cycles(idx1:idx2) = c;
    
                  % Assign difficulty to this interval
                  if params.sim.epochs.difficulty.apply
                      overtake_info.difficulty(idx1:idx2) = get_difficulty(seq_diff, data.sim.sim2track.matrix, ...
                                                               mean([this_left this_right]), ...
                                                               params.sim.lane_dist);
                  end
                  
                  % Assign outcome to this interval
                  if params.sim.epochs.outcomes.apply
                      outcome = get_outcome(data.sim.sim2track, ...
                                                               [this_left this_right]);
                      if outcome > 0
                         outcome = 1;
                      elseif outcome < 0
                         outcome = -1;
                      end
                      overtake_info.outcomes(idx1:idx2) = repmat(outcome,idx2-idx1+1,1);
                  end
              
              end
                                                   
          end
          % Reset
          this_left = -1;
       end
    end

end
