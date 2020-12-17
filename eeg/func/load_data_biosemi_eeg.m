function [ data ] = load_data_biosemi_eeg ( params, subject )
% 
% Loads data from a subject's EEG files from Biosemi format (bdf). 
% Finds simulation start trigger and zeros time series on this data
% point.
%

% Check that data exists; unzip if necessary
data = [];

subj_dir = sprintf('%s/%s/%s', params.io.input_dir, params.eeg.sub_dir, subject);

cfg = params.eeg.cfg;
cfg.headerfile = sprintf('%s/%s%s%s.bdf', subj_dir, params.eeg.convert.prefix, subject, params.eeg.convert.suffix);
if ~exist(cfg.headerfile, 'file')
   warning('Subject %s has no EEG data. Skipping...', subject);
   return;
end

cfg.dataset = cfg.headerfile;
[~,data.eeg.ft] = evalc('ft_preprocessing(cfg)');

data.eeg.eeg_channels = [];
data.eeg.eog_channels = [];
data.eeg.all_channels = data.eeg.ft.label;
data.eeg.idx_eog = [];
for c = 1 : length(data.eeg.ft.label)
   channel = data.eeg.ft.label{c};
   if contains(channel,'EOG')
       data.eeg.eog_channels = [data.eeg.eog_channels {channel}];
       data.eeg.idx_eog = [data.eeg.idx_eog c];
   else
       data.eeg.eeg_channels = [data.eeg.eeg_channels {channel}];
   end
end

marker_file = cfg.headerfile;
[~,cfg.event] = evalc('ft_read_event(marker_file)');

start_idx =[];
restart_time = 0;

if ~isempty(params.eeg.convert.restart_file)
    input_file = sprintf('%s/%s/%s', params.io.input_dir, ...
                                     params.io.metadata_dir, ...
                                     params.eeg.convert.restart_file);
                         
    opts = detectImportOptions(input_file);
    opts.VariableTypes(1)={'char'};
    opts.VariableNames(1)={'SubjectID'};
    T_restart = readtable(input_file, opts);

    idx = find(strcmp(T_restart.SubjectID, subject),1);
    if ~isempty(idx)
        restart_time = T_restart.RestartTime(idx);
        warning(' %s: Restarting at %1.2fs.', subject, restart_time);
        [~,start_idx] = min(abs(data.eeg.ft.time{1}-restart_time));
        triggers = get_struct_values( cfg.event, 'value', true );
        start_idx = find(triggers >= start_idx,1);
    end
end

if isempty(start_idx) && ~isempty(params.eeg.convert.start_byte)
    % Look for start value
    triggers = get_struct_values( cfg.event, 'value', true );
    start_idx = find(triggers == params.eeg.convert.start_byte, 1); 
    if ~isempty(start_idx)
        restart_time = data.eeg.ft.time{1}(cfg.event(start_idx).sample);
    end
end

if ~isempty(start_idx)
    cfg2 = [];
    cfg2.latency = [restart_time data.eeg.ft.time{1}(end)];
    [~,data.eeg.ft] = evalc('ft_selectdata(cfg2, data.eeg.ft);');
    data.eeg.ft.time(1) = {data.eeg.ft.time{1} - restart_time};

    tidx = get_struct_values(cfg.event, 'sample', true);
    for i = 1 : length(cfg.event)
        cfg.event(i).sample = cfg.event(i).sample - tidx(start_idx) + 1;
    end
    cfg.event = cfg.event(start_idx:end); 
end

% -- Process EEG events, adjust time to be relative to first event --
is_stim = strcmp({cfg.event(:).type}, 'STATUS');

tidx = get_struct_values(cfg.event, 'sample', true);
tidx = tidx(is_stim);
value = get_struct_values(cfg.event, 'value', true);
value = value(is_stim);

triggers = value;
type = get_struct_values(cfg.event, 'type');

if params.eeg.convert.align_sim_triggers
    % Align triggers in EEG to simulation log times
    % Use this only if the EEG trigger values are incorrect
    event_file = sprintf('%s/%s/sim/events-All.csv', params.io.output_dir, subject);
    opts = detectImportOptions(event_file);
    opts.VariableTypes = {'double', 'double', 'double', 'char', 'double'};
    T_events_sim = readtable(event_file);
    
    event_types = T_events_sim.EventType;
    idx_start = find(strcmp(event_types, 'SimulatorStarted'),1);
    if ~isempty(idx_start) && idx_start > 0
        % Remove any events before start of simulation
        T_events_sim = T_events_sim(idx_start:end,:);
    end
    
    byte_sim = T_events_sim.AdjSerialByte;
    idx_rm = isnan(byte_sim);
    t_sim = T_events_sim.Millis / 1000;
    t_sim = t_sim - t_sim(1);
    t_sim = t_sim(~idx_rm);
    byte_sim = byte_sim(~idx_rm);
    
    t_eeg = data.eeg.ft.time{1};
    t_eeg = t_eeg(tidx);
    t_eeg = t_eeg - t_eeg(1);
    
    [idx_map,d] = align_timeseries_eeg(t_eeg, t_sim);
    triggers = byte_sim(idx_map(:,2));
    
    type = type(idx_map(:,1));
    tidx = tidx(idx_map(:,1));
    value = value(idx_map(:,1));
    
else
    % Add 256 to each reset of byte values (SerialByte = 1)
    % This is because triggers increment at each step, and reset at 
    % 255 since they are 8-bit values.
    to_add = 0;
    first_found = false;
    for i = 1 : length(triggers)
        if triggers(i) == 1
            if first_found
                to_add = to_add + 256;
            else
                first_found = true;
            end
        end
        triggers(i) = triggers(i) + to_add;
    end
    
end


times = data.eeg.ft.time{1}(tidx);
data.eeg.trigger_idx = tidx;
data.eeg.cfg = cfg;
data.eeg.events = table(times(:), tidx(:), triggers(:), value(:), type(:));
data.eeg.events.Properties.VariableNames = {'Time', 'Index', 'Trigger', 'Value', 'Type'};



    

end