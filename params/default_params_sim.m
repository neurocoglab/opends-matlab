%% This script defines default simulation log parameters for the opends 
% processing pipeline

%% General
params.sim.sub_dir = 'driving';
params.sim.lane_dist = 7035.461;
params.sim.version = '1.8.4';

%% Convert log files
params.sim.convert.prefix = 'simlog-';
params.sim.convert.filter = 'serial-filter';
params.sim.convert.exec = 'simlog2csv';


%% Pre-process rounds
params.sim.rounds.plots.save = true;
params.sim.rounds.fix_repeats = true;
params.sim.rounds.plots.color = [0.8 0.8 0.9];

params.sim.log.trigger_format = '%d64 %f %d %q %q %q %q %f %f %d %d %d %d';
params.sim.log.lanechange_format = '%d64 %f %d %q %q %q %q %q %d %f %d %d %f %d %d';
params.sim.log.accelerate_format = '%d64 %f %d %q %q %f %f %f %d %d %d %f %d %d';
params.sim.log.reward_format = '%d64 %f %d %q %q %q %f %d %f %d %d';
params.sim.log.state_format = '%d64 %f %d %q %d %d %f %d %f %d %d';
params.sim.log.start_format = '%d64 %f %d %q %f %d %d';
params.sim.log.events_format = '%d64 %f %d %q %d';
params.sim.log.resumed_format = '%d64 %f %d %q %q %d %d %d %d %d %f';
params.sim.log.simended_format = '%d64 %f %d %q %f %d %d';
params.sim.log.rewarddisp_format = '%d64 %f %d %q %q %q %d %d %f %d %d';
params.sim.log.messagebutton_format = '%d64 %f %d %q %q %q %d %q %f %d %d %d %d';
params.sim.log.roadsign_format = '%d64 %f %d %q %q %q %f %d %d';

% Fix for distances (bug in early versions)
if exist('fixdist.mat', 'file')
    params.sim.rounds.fixdist = load('fixdist.mat');
end

% String pattern for max cycles (implicit iterate-cycle)
params.sim.log.max_cycle_str = 'Max cycles reached.'; % TODO: make this max-cycles-reached

params.sim.baseline.intervals_file = 'baseline_intervals.csv';

% Message button presses
params.sim.rounds.messagebutton.apply = false;
params.sim.rounds.messagebutton.maps_file = 'message_button_map.csv';

% Roadsign changes
params.sim.rounds.roadsign.apply = false;
params.sim.rounds.roadsign.maps_file = 'roadsign_map.csv';

%% Metadata
params.sim.metadata.uid_field = 'UID';
params.sim.metadata.score_field = 'Game_score';

%% Process epochs
params.sim.epochs.difficulty.apply = true;
params.sim.epochs.difficulty.sequence_file = 'sequence_difficulty.csv';
params.sim.epochs.difficulty.levels = [1 2];
params.sim.epochs.difficulty.labels = {'Easy','Difficult'};
params.sim.epochs.outcomes.apply = true;
params.sim.epochs.outcomes.levels = [-1 1];

%% Process events
params.sim.events.difficulty.apply = true;
params.sim.events.difficulty.levels = [1 2];
params.sim.events.outcomes.apply = true;
params.sim.events.outcomes.levels = [-1 1];
params.sim.events.difficulty.labels = {'Easy','Difficult'};
params.sim.events.button_presses.apply = true;

% Includes traffic decision button presses?
params.sim.events.traffic_decision.apply = true;
params.sim.events.traffic_decision.write_to_file = true;

% Includes road sign changes?
params.sim.events.roadsign.apply = true;

% Includes fixations to areas of interest?
params.sim.events.fixations.apply = false;
params.sim.events.fixations.aoi_file = '';

%% Plots
params.sim.plots.pass_color = [0.9 0.9 0.95];
params.sim.plots.baseline_color = [0.9 0.95 0.9];



