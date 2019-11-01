%% Read EEG for driving study

eeg_dir = '/Volumes/AndrewElements/data/driving/eeg';
subject = '4482';

addpath '/Users/lpzatr/Documents/MATLAB/lib/fieldtrip-20180110'

%% Read in sim log and match EEG markers with log markers

% To temporally align all time series, express time for all observations 
% (tracker, EEG, simlog) relative to simulation start event

cfg = [];
cfg.headerfile = sprintf('%s/%s/%s-eeg/%s.vhdr', eeg_dir, subject, subject, subject);
cfg.dataset = sprintf('%s/%s/%s-eeg/%s.eeg', eeg_dir, subject, subject, subject);
cfg.continuous = true;
cfg.blocksize = 60;

[data] = ft_preprocessing(cfg);

cfg.event = ft_read_event(sprintf('%s/%s/%s-eeg/%s.vmrk', eeg_dir, subject, subject, subject));

preproc = load('preproc_params_hd.mat');

% Match EEG times to sim log events

is_stim = strcmp({cfg.event(:).type}, 'Stimulus');
tidx = cell2mat({cfg.event(is_stim).sample});
value = {cfg.event(is_stim).value};
type = {cfg.event(is_stim).type};

triggers = zeros(1,length(value));

for i = 1 : length(value)
    event = value{i};
    triggers(i) = str2num(event(2:end));
end

% First trigger is simulation start
start_sim = tidx(1);
t_simstart = data.time{1}(start_sim);
triggers = triggers(2:end);
value = value(2:end);
type = type(2:end);

% Add 256 to each reset of byte values (SerialByte = 1)
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

times = data.time{1}(tidx);
times = times - t_simstart; % times(1);
times = times(2:end);
tidx = tidx(2:end);

eeg.events = table(times', triggers', value', type');
eeg.events.Properties.VariableNames = {'Time', 'Trigger', 'Value', 'Type'};

outdir = sprintf('%s/%s/%s', preproc.params.root_dir, preproc.params.output_dir, subject);

% Get timestamp of simulation start
input_file = sprintf('%s/events-SimulatorStarted.csv', outdir);
[values, hdr] = import_log(input_file, preproc.params.simlog.start_format);
log_start_ms = values{1};

% Load lange change events
input_file = sprintf('%s/events-LaneChange.csv', outdir);
[values, hdr] = import_log(input_file, preproc.params.simlog.lanechange_format);

T = values{1};
mylog.times = double(T-log_start_ms) / 1000;
mylog.triggers = values{find(strcmp(hdr,'AdjSerialByte'))};
mylog.logid = values{find(strcmp(hdr,'LogId'))};
mylog.type = values{find(strcmp(hdr,'EventType'))};
mylog.lane_from = values{find(strcmp(hdr,'LaneFrom'))};
mylog.lane_to = values{find(strcmp(hdr,'LaneTo'))};
mylog.sim_time =values{find(strcmp(hdr,'SimulationTime'))};
mylog.cycle=values{find(strcmp(hdr,'Sim:Game:Cycle'))};
mylog.repeat=values{find(strcmp(hdr,'Sim:Game:Repeat'))};

sim.events = table(mylog.times, mylog.triggers, mylog.logid, mylog.type, mylog.lane_from, ...
                   mylog.lane_to, mylog.sim_time, mylog.cycle, mylog.repeat);
sim.events.Properties.VariableNames = {'Time', 'Trigger', 'LogId', 'Type', 'LaneFrom', ...
                                       'LaneTo', 'SimTime', 'Cycle', 'Repeat'};

writetable(eeg.events, sprintf('eeg-events-%s.csv', subject));
writetable(sim.events, sprintf('sim-events-%s.csv', subject));

clear mylog;
clear data;
clear eeg;
clear sim;

%% Eye tracker events

preproc = load('preproc_params_hd.mat');

outdir = sprintf('%s/%s/%s', preproc.params.root_dir, preproc.params.output_dir, subject);
results_file = sprintf('%s/results.mat',outdir);
flag_file = sprintf('%s/sim_logs.done', outdir);

if ~exist(flag_file, 'file')
   error('ET logs have not been preprocessed for subject "%s".\n', subject)
end

load(results_file);

% Get first event time [tracker time AND log time]; get sim start time [log time]
% Compute sim start time as tracker time; subtract from time series

logid = 1;

t_simstart_sim = data.sim.events.values{find(strcmp(data.sim.events.hdr,'Millis'))}(1);
ids = data.sim.events.values{find(strcmp(data.sim.events.hdr,'LogId'))};
t_event_sim = data.sim.events.values{find(strcmp(data.sim.events.hdr,'Millis'))}(find(ids,logid));
t_delta = t_event_sim - t_simstart_sim; % Time between sim start and first event, in ms
ids = data.eye.log.messages{find(strcmp(data.eye.log.hdr,'LogId'))};
t_event_tracker =  data.eye.log.messages{find(strcmp(data.eye.log.hdr,'Time'))}(find(ids,logid));

t_event_tracker = t_event_tracker / 1000; % tracker unit is microseconds, convert to ms
t_simstart_tracker = t_event_tracker - t_delta; % simulation start in tracker time

t_delta = single(data.eye.t_start - t_simstart_tracker); % Time is currently relative to first timepoint,
                                                         % make it relative to sim start 

results.t = results.t + t_delta;
results.t = results.t / 1000;  % Convert to seconds

% Reload EEG
cfg = [];
cfg.headerfile = sprintf('%s/%s/%s-eeg/%s.vhdr', eeg_dir, subject, subject, subject);
cfg.dataset = sprintf('%s/%s/%s-eeg/%s.eeg', eeg_dir, subject, subject, subject);
cfg.continuous = true;
cfg.blocksize = 60;
cfg.hpfilter = 'no';
cfg.hpfreq = 1;
cfg.lpfilter = 'no';
cfg.lpfreq = 100;

[data.eeg] = ft_preprocessing(cfg);
data.eeg.time = {data.eeg.time{1} - t_simstart};


blah = 0; % Prevents editor from being stupid while editing

%% Plot stuff

h=figure;

plot(results.t/60,zscore(results.diam_left)+5,'g');
hold on;
plot(results.t/60,zscore(results.pos_left_x)/2+2,'r');
plot(results.t/60,zscore(results.pos_left_y)/2-1,'b');

channel = 'hEOGleft';
plot(data.eeg.time{1}/60, zscore(data.eeg.trial{1}(find(strcmp(data.eeg.label,channel)),:))/2-6, 'k');

xlim([0 35]);
ylim([-10,8]);

blah = 0;



%% 

proc = load('processing_params.mat');
process_results_file = sprintf('%s/%s/processing_results.mat',proc.params.data_dir,subject);

proc.results = load(process_results_file);


