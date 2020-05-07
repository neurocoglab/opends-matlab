function [ trials ] = get_trials_eeg( data, sim, params )

% Returns a set of trials, and associated baseline trials, dervied from an 
% EEG dataset.
%
% params should be a sub-struct with "baseline" and "windows" as fields (e.g.,
% pass "params.eeg.erp").
%
% Currently only processes LaneChangeLeft and LaneChangeRight events
%


trials = [];   
Ts = 1000 / data.eeg.ft.fsample; % Sample period in ms
baseline_offset = params.baseline.offset; % ms
baseline_offset = round(baseline_offset / Ts); % samples
time_tol = 0.01;
baseline_max = params.baseline.maxlen;
baseline_min = params.baseline.minlen;

bw_window = params.windows('Baseline',:);

% LaneChangeLeft events
T = sim.events(strcmp(sim.events.LaneFrom,'Lane.hiway.1'),{'Time','Trigger','Difficulty','Outcomes','Cycle','Repeat'});
window = params.windows('LaneChangeLeft',:);
[trials.left_change.trl, trials.left_change.trl_baseline] = get_trials(T, window, bw_window);

% LaneChangeRight events
T2 = sim.events(strcmp(sim.events.LaneFrom,'Lane.hiway.2'),{'Time','Trigger','Difficulty','Outcomes','Cycle','Repeat'});
window = params.windows('LaneChangeRight',:);
[trials.right_change.trl, trials.right_change.trl_baseline] = get_trials(T2, window, bw_window, T);


    function [trl, trl_baseline] = get_trials(T, window, bw_window, T_bl)
        
        if nargin < 4
            T_bl=[];
        end
        
        trl = zeros(0,8);
        trl_baseline = zeros(0,8);
        trial_idx = window.TrialIdx;
        window = [window.From / Ts, window.To / Ts];
        bw_window = [bw_window.From / Ts, bw_window.To / Ts];
        
        skipped = 0;

        % Define trials
        for i = 1 : height(T)
            record_i = data.eeg.events(data.eeg.events.Trigger==T.Trigger(i),:);
            if isempty(record_i)
                warning('   No EEG event matching trigger %d...', T.Trigger(i));
            else
                ok=1;
                % Deal with double triggers (should not occur but do)
                if height(record_i) > 1
                    t_eeg = record_i{:,{'Time'}};
                    t_sim = T{i,{'Time'}};
                    ok = abs(t_eeg - t_sim) < time_tol;
                end
                for j = 1 : length(ok)
                    if ok(j)
                        record_ij = record_i(j,:);
                        window_ij = [record_ij.Index + window(1), record_ij.Index + window(2)];
                        bw_window_ij = [];
                        
                        if ~isempty(T_bl)
                            % Offset from pre-defined baseline trial window
                           dt = T.Time(i) - T_bl.Time;
                           idx = find(dt>0);
                           if ~isempty(idx) && dt(idx(end)) <= baseline_max && dt(idx(end)) >= baseline_min
                               bw_record_ij = data.eeg.events(data.eeg.events.Trigger==T_bl.Trigger(idx(end)),:);
                               bw_window_ij = [bw_record_ij.Index + bw_window(1), bw_record_ij.Index + bw_window(2)];
                               bw_window_ij = bw_window_ij + baseline_offset;
                           else
                               if isempty(idx)
                                   fprintf('DEBUG: no baseline found for trial.\n');
                               else
                                   fprintf('DEBUG: no start event found for trial (dt=%1.3f).\n', dt(idx(end)));
                               end
                           end
                        else
                            % Offset from current trial window
                            bw_window_ij = [record_ij.Index + bw_window(1), record_ij.Index + bw_window(2)] + baseline_offset;
                        end
                        
                        if isempty(bw_window_ij) || bw_window_ij(1) + min(window(1),0) < 0 || window_ij(2) > length(data.eeg.ft.time{1})
                            skipped = skipped + 1;
                        else
                            trl(end+1,:) = [window_ij(1), window_ij(2), window(1), trial_idx, T.Difficulty(i), ...
                                            T.Outcomes(i), T.Cycle(i), T.Repeat(i)];

                            if bw_window_ij(1) > 0 && bw_window_ij(2) < length(data.eeg.ft.time{1})
                                trl_baseline(end+1,:) = [bw_window_ij(1), bw_window_ij(2), window(1), ...
                                                         trial_idx, T.Difficulty(i), T.Outcomes(i), T.Cycle(i), T.Repeat(i)];
                            else
                                trl_baseline(end+1,:) = nan(6,1);
                            end
                        end
                    end
                end
            end
        end

        trl(isnan(trl)) = 0;
        trl_baseline(isnan(trl_baseline)) = 0;
        
        if skipped > 0
           fprintf(' ** Skipped %d trials whose baselines could not be determined, or which extended beyond end of time series.\n', skipped);
        end
        
    end

end