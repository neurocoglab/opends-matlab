function [ data ] = load_data_brainvision_eeg ( params, subject )
% 
% Loads data from a subject's EEG files from Brain Vision format (vhdr and
% vmrk). Finds simulation start trigger and zeros time series on this data
% point.
%

ok = true;

% Check that data exists; unzip if necessary
data = [];

subj_dir = sprintf('%s/%s/%s', params.io.input_dir, params.eeg.sub_dir, subject);

cfg = params.eeg.cfg;
cfg.headerfile = sprintf('%s/%s-eeg/%s.vhdr', subj_dir, subject, subject);
if ~exist(cfg.headerfile, 'file')
   zip_file = sprintf('%s/%s-eeg.zip', subj_dir, subject);
   unzip_dir = sprintf('%s/%s-eeg', subj_dir, subject);
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

cfg.dataset = sprintf('%s/%s-eeg/%s.eeg', subj_dir, subject, subject);
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

marker_file = sprintf('%s/%s-eeg/%s.vmrk', subj_dir, subject, subject);
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


end