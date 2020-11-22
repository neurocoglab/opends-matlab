function [ trials, baseline_trials ] = get_sim_trials_eeg( params, data )
%   GET_SIM_TRIALS_EEG Produce trial definitions based on simulation
%   events and triggers in the EEG data
%      

trials = {};
baseline_trials = {};

W = readtable(params.eeg.trials.windows_file);

% LeftChange events
idx = find(strcmp(W.Event,'ChangeLeft'));

if ~isempty(idx)
    T = data.sim.lane_change.values(strcmp(data.sim.lane_change.values.LaneFrom,'Lane.hiway.1'),:);
    window = W(idx,:);
    [trl,trl_baseline] = get_trials(data, T, window);
    trials(end+1) = {trl};
    baseline_trials(end+1) = {trl_baseline};
end

% RightChange events
T = sim.events(strcmp(sim.events.LaneFrom,'Lane.hiway.2'),{'Time','Trigger'});
window = W('LaneChangeRight',:);
[trl,trl_baseline] = get_trials(T, window);
trials(end+1) = {trl};
baseline_trials(end+1) = {trl_baseline};

% ButtonPress events


function [trl,trl_baseline] = get_trials(data, T, window)

    trial_idx = window.TrialIdx;

    % For every event in T, add a trial with a window around it
    % Also add index, difficulty, and outcome
    
    
    
    % Overtake difficulty
    process_results_file = sprintf('%s/%s/results_eye.mat', params.io.output_dir, data.subject);
    processing_results = load(process_results_file);
    diffs = processing_results.results.events.left_change.diffs;

    Ts = 1 / data.eeg.ft.fsample * 1000; % Sample period in ms
    window = [window.From / Ts, window.To / Ts];
    trl = zeros(0,5);
    trl_baseline = zeros(0,5);
    baseline_offset = -20.0; % seconds
    baseline_offset = round(baseline_offset * 1000 / Ts);

    time_tol = 0.01;

    % Define Left-lane-change trials
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
                    trl(end+1,:) = [window_ij(1), window_ij(2), window(1), trial_idx, diffs(i)]; 
                    bw_window = window_ij + baseline_offset;
                    if bw_window(1) > 0 && bw_window(2) < length(data.eeg.ft.time{1})
                        trl_baseline(end+1,:) = [bw_window(1), bw_window(2), window(1), trial_idx, diffs(i)];
                    else
                        trl_baseline(end+1,:) = nan(5,1);
                    end
                end
            end
        end
    end

    idx_nan = find(any(isnan(trl_baseline),2));
    if ~isempty(idx_nan)
       trl(idx_nan,:) = [];
       trl_baseline(idx_nan,:) = [];
    end


end




end

