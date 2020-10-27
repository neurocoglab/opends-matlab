function [ result ] = get_timefreq_trials_eeg( params, data, results )
%get_timefreq_trials Produces trials from rounds (and repeats),
%   suitable for time/frequency analysis. 
%
% N_rounds = # rounds + # repeats
% N_types  = # event types
%
% result.time:         1 x N_rounds cell array of sample times, in seconds,
%                       zeroed at the start of the round, with each cell
%                       containing an N_timepoints array
% result.trial:        1 x N_rounds cell array of EEG signals, with
%                       each cell containing an N_channels X N_timepoints
%                       matrix
% result.event:       N_types x N_rounds cell array of event indices, with
%                       each cell containing a list of indices
% result.event_table:  Table of event types and phases
%


% Get baseline relative to simulation start
baseline_times = data.sim.sim2track.baseline / 1000;
baseline_times = mean(baseline_times,2);

% Get EEG trigger for rounds
T_rounds = data.sim.rounds.triggers;
T_start = T_rounds(strcmp(T_rounds.Position, 'start'),:);
T_stop = T_rounds(strcmp(T_rounds.Position, 'stop'),:);

result = [];
result.fsample = data.eeg.ft.fsample;
result.hdr = data.eeg.ft.hdr;

T_event_types = readtable(params.eeg.timefreq.event_types_file);
event_types = unique(T_event_types.EventType);
result.event_table = T_event_types;

N_event_info = 0;
result.event_labels = {};
for j = 1 : length(event_types)
    idx_et = find(strcmp(T_event_types.EventType, event_types{j}));
    event_info = unique(T_event_types.Phase(idx_et));
    N_event_info = N_event_info + length(event_info);
    for kk = 1 : length(event_info)
        result.event_labels = [result.event_labels;{event_types{j}, event_info{kk}}];
    end
end

N_rounds = height(T_start);
N_types = length(event_types);
result.time = cell(1, N_rounds);
result.baseline = cell(1, N_rounds);
result.trial = cell(1, N_rounds);
result.event = cell(N_types, N_rounds);
result.event_info = cell(N_event_info, N_rounds);

result.label = data.eeg.ft.label;
result.rounds = zeros(1,N_rounds);
result.repeats = zeros(1,N_rounds);

buffer = params.eeg.timefreq.buffer * data.eeg.ft.fsample;

% For each round/repeat
for i = 1 : height(T_start)
    %  Get start/end indices +/- buffer
    trigger_start = T_start.AdjSerialByte{i};
    trigger_end = T_stop.AdjSerialByte{i};
    
    %  Get corresponding times, zero @ start_index+buffer
    idx_start = find(data.eeg.events.Trigger == trigger_start, 1);
    idx_end = find(data.eeg.events.Trigger == trigger_end, 1);
    isok = true;
    if ~isempty(idx_start)
        idx_start = data.eeg.events.Index(idx_start);
    else
        % Deal with missing trigger
        if params.general.debug 
            fprintf('\nDEBUG: EEG missing start trigger %d\n', trigger_start);
        end
        idx_start = get_trigger_index( data.sim.events.values, trigger_start );
        isok = ~isempty(idx_start);
    end
    if ~isempty(idx_end)
        idx_end = data.eeg.events.Index(idx_end);
    else
        % Deal with missing trigger
        if params.general.debug 
            fprintf('\nDEBUG: EEG missing end trigger %d\n', trigger_end);
        end
        idx_end = get_trigger_index( data.sim.events.values, trigger_end );
        isok = ~isempty(idx_start);
    end
    
    if isok
        
        % Get baseline periods in this window
        time = data.eeg.ft.time{1}(idx_start:idx_end);
        idx_bl = find(baseline_times > time(1) & baseline_times < time(end));
        result.baseline(i) = {baseline_times(idx_bl) - data.eeg.ft.time{1}(idx_start)};
    
        %  Slice corresponding data
        idx1 = max(1, idx_start-buffer);
        idx2 = min(length(data.eeg.ft.time{1}), idx_end+buffer);
        
        time = data.eeg.ft.time{1}(idx1:idx2);
        time = time - data.eeg.ft.time{1}(idx_start);
        result.time(i) = {time};
        result.trial(i) = {data.eeg.ft.trial{1}(:,idx1:idx2)};
        result.rounds(i) = T_start.Cycle{i};
        result.repeats(i) = T_start.Repeat{i};
        
        %  For each event, get indices relative to round start
        T_events = data.eeg.events(data.eeg.events.Index >= idx_start & ...
                                   data.eeg.events.Index <= idx_end, :);
        T_events.Index = T_events.Index - idx_start + 1;
        
        idx_p = 0;
        
        for j = 1 : length(event_types)
            event_type = event_types{j};
            idx_et = find(strcmp(T_event_types.EventType, event_type));
            event_info = unique(T_event_types.Phase(idx_et));
            [events, phases] = get_events( event_type, T_events, event_info );
            idx_p = idx_p(end)+1:idx_p(end)+length(phases);
            result.event(j,i) = {time(events)};
            result.event_info(idx_p,i) = phases;

            % Binarize pos/neg phases
            for l = 1 : length(event_info)
                T_p = T_event_types(strcmp(T_event_types.EventType, event_type) & ...
                                    strcmp(T_event_types.Phase, event_info{l}),:);
                levels = T_p.Level;
                if any(strcmp(levels, 'negative')) && any(strcmp(levels, 'positive'))
                    phases = result.event_info{l,i};
                    P = phases;
                    P(phases < 0) = T_p.Index(strcmp(levels, 'negative'));
                    P(phases > 0) = T_p.Index(strcmp(levels, 'positive'));
                    result.event_info(idx_p(l),i) = {P};
                end
            end
        end
    
    end
end

result.eeg.timefreq = result;


    % Get index in EEG time series of closest previous trigger in the
    % EEG trigger list; this compensates for missing triggers
    function idx = get_trigger_index( T_events, target_trigger )
        
        time_t = T_events.Millis(find(T_events.AdjSerialByte==target_trigger,1));
        tr_prev = target_trigger - 1;
        idx = [];
        
        while tr_prev > 0
            ii = find(data.eeg.events.Trigger == tr_prev, 1);
            if ~isempty(ii)
                time_prev = T_events.Millis(ii);
                delta_t = (time_t - time_prev) / 1000;
                time_eeg = double(data.eeg.events.Time(ii) + delta_t);
                % Get closest index
                [ ~, idx ] = min( abs( double(data.eeg.ft.time{1}) - time_eeg ) );
                return;
            end
            tr_prev = tr_prev - 1;
        end
        
    end


    function [events, phases] = get_events( event_type, T_events, event_info )
       
        events = [];
        [T_sim_events, phases] = get_sim_events( event_type, event_info );
        
        if isempty(T_sim_events)
           return; 
        end
        
        idx = [];
        for k = 1 : height(T_events)
            trigger = T_events.Trigger(k);
            idx_k = find(T_sim_events.AdjSerialByte == trigger,1);
            if ~isempty(idx_k)
                events = [events T_events.Index(k)];
                idx = [idx idx_k];
            end
        end
        
        for k = 1 : length(phases)
            phases(k) = {phases{k}(idx)};
        end
    end

    % Get all the simulation events for the given type,
    % and populate the "phases" for each event as well
    % "event_type" and "event_info" must match the types
    % produced by "run_processing_eye_sim.m"; i.e.,
    % in "results.eye.events.(event_type).(phase_type)"
    function [table, phases] = get_sim_events( event_type, event_info )
        
        table = [];
        phases = [];
        
        if strcmp(event_type, 'left_change') == 1
            table = data.sim.lane_change.values;
            table = table(find(strcmp(table.Direction, 'left')==1),:);
        elseif strcmp(event_type, 'right_change') == 1
            table = data.sim.lane_change.values;
            table = table(strcmp(table.Direction, 'right')==1,:);
        else
            warning('Invalid event type: %s', event_type);
            return;
        end
        
        phases = cell(1,length(event_info));
        for k = 1 : length(event_info)
           phases(k) = {results.eye.events.(event_type).(event_info{k})}; 
        end
        
    end

end

