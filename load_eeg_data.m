function [ data ] = load_eeg_data ( params, preproc, subject )

ok = true;

% Find bad channels for this subject, if defined
bad_channels = [];

if exist(params.eeg.bad_channel_file, 'file')
   opts = detectImportOptions(params.eeg.bad_channel_file);
   opts = setvartype(opts, 'char');
   bad_channels = readtable(params.eeg.bad_channel_file, opts);
   bad_channels.Properties.VariableNames = [{'Subject'},{'Channel'}];
end

% Check that data exists; unzip if necessary
data = [];

cfg = params.eeg.cfg;
cfg.headerfile = sprintf('%s/%s/%s-eeg/%s.vhdr', params.eeg.data_dir, subject, subject, subject);
if ~exist(cfg.headerfile, 'file')
   zip_file = sprintf('%s/%s/%s-eeg.zip', params.eeg.data_dir, subject, subject);
   unzip_dir = sprintf('%s/%s/%s-eeg', params.eeg.data_dir, subject, subject);
%        fprintf('Zip file: %s\n', zip_file);
   if ~exist(zip_file, 'file')
       warning('Subject %s has no EEG data. Skipping...', subject);
       ok = false;
   else
       unzip(zip_file, unzip_dir);
       ok = exist(cfg.headerfile, 'file');
   end
end

% Load data
if ~ok, return; end

cfg.dataset = sprintf('%s/%s/%s-eeg/%s.eeg', params.eeg.data_dir, subject, subject, subject);
[~,data.eeg.ft] = evalc('ft_preprocessing(cfg)');

data.eeg.badchannels = [];

% Remove bad channels
if ~isempty(bad_channels)
    T = bad_channels{strcmp(bad_channels.Subject, subject),2};
    data.eeg.badchannels = T;
    if ~isempty(T)
        cfg2 = [];
        cfg2.channel = data.eeg.ft.label;
        idx_rem = [];
        for c = 1 : length(T)
            idx_rem = [idx_rem find(strcmp(cfg2.channel,T{c}))];
%             cfg2.channel(find(strcmp(cfg2.channel,T{c})))=[];
        end
        cfg2.channel(idx_rem) = [];
        [~,data.eeg.ft] = evalc('ft_selectdata(cfg2, data.eeg.ft)');
        
        fprintf('%s: Removed %d bad channels.\n', subject, length(T));
    end
end

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

marker_file = sprintf('%s/%s/%s-eeg/%s.vmrk', params.eeg.data_dir, subject, subject, subject);
[~,cfg.event] = evalc('ft_read_event(marker_file)');
idx_eeg_start = find(strcmp(vertcat({cfg.event.value}),'S128'));
idx_eeg_start = idx_eeg_start(1);
data.eeg.t_eeg_start = data.eeg.ft.time{1}(idx_eeg_start);

% -- Process EEG events, adjust time to be relative to first event --

is_stim = strcmp({cfg.event(:).type}, 'Stimulus');
tidx = cell2mat({cfg.event(is_stim).sample});
value = {cfg.event(is_stim).value};
type = {cfg.event(is_stim).type};
triggers = zeros(1,length(value));

for i = 1 : length(value)
    event = value{i};
    triggers(i) = str2num(event(2:end));
end

% First trigger is simulation start; set times relative to this
start_sim = tidx(1);
t_simstart_eeg = data.eeg.ft.time{1}(start_sim);
triggers = triggers(2:end);
value = value(2:end);
type = type(2:end);

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

times = data.eeg.ft.time{1}(tidx);
times = times - t_simstart_eeg;
times = times(2:end);
tidx = tidx(2:end);

data.eeg.t_simstart = t_simstart_eeg;
data.eeg.trigger_idx = tidx;
data.eeg.cfg = cfg;
data.eeg.events = table(times', tidx', triggers', value', type');
data.eeg.events.Properties.VariableNames = {'Time', 'Index', 'Trigger', 'Value', 'Type'};

data.eeg.ft.time = {data.eeg.ft.time{1} - t_simstart_eeg};

% Get simulation time start for later synching
outdir = sprintf('%s/%s/%s', preproc.params.root_dir, preproc.params.output_dir, subject);
results_file = sprintf('%s/results.mat',outdir);
flag_file = sprintf('%s/sim_logs.done', outdir);

if ~exist(flag_file, 'file')
   error('Simulation logs have not been preprocessed for subject "%s".\n', subject)
end

results_eye = load(results_file);
logid = 1;

% t_simstart_sim = data.sim.simstart.values{find(strcmp(data.sim.events.hdr,'Millis'))}(1);
% ids = data.sim.events.values{find(strcmp(data.sim.events.hdr,'LogId'))};
% t_event_sim = data.sim.events.values{find(strcmp(data.sim.events.hdr,'Millis'))}(find(ids==logid));
% t_delta = t_event_sim - t_simstart_sim; % Time between sim start and first event, in ms
% ids = data.eye.log.messages{find(strcmp(data.eye.log.hdr,'LogId'))};
% t_event_tracker =  data.eye.log.messages{find(strcmp(data.eye.log.hdr,'Time'))}(find(ids==logid));

t_simstart_sim = results_eye.data.sim.simstart.values{find(strcmp(results_eye.data.sim.events.hdr,'Millis'))}(1);
ids = results_eye.data.sim.events.values{find(strcmp(results_eye.data.sim.events.hdr,'LogId'))};
t_event_sim = results_eye.data.sim.events.values{find(strcmp(results_eye.data.sim.events.hdr,'Millis'))}(find(ids==logid));
t_delta = t_event_sim - t_simstart_sim; % Time between sim start and first event, in ms


data.eeg.t_delta_sim = t_delta;
data.eeg.t_start_sim = t_simstart_sim;


end