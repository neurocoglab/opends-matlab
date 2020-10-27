%% This script defines default simulation log parameters for the opends 
% processing pipeline

%% General
params.sim.sub_dir = 'driving';
params.sim.lane_dist = 7035.461;

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
params.sim.log.resumed_format = '%d64 %f %d %q %q %d %d %d %d %d';
params.sim.log.simended_format = '%d64 %f %d %q %f %d %d';
params.sim.log.rewarddisp_format = '%d64 %f %d %q %q %q %d %d %f %d %d';
params.sim.log.messagebutton_format = '%d64 %f %d %q %q %q %d %q %f %d %d %d %d';
params.sim.log.roadsign_format = '%d64 %f %d %q %q %q %f %d %d';

% Fix for distances (bug in early versions)
if exist('fixdist.mat', 'file')
    params.sim.rounds.fixdist = load('fixdist.mat');
end

params.sim.baseline.intervals_file = 'baseline_intervals.csv';

% Message button presses
params.sim.rounds.messagebutton.apply = false;
params.sim.rounds.messagebutton.maps_file = 'message_button_map.csv';

% Roadsign changes
params.sim.rounds.roadsign.apply = false;
params.sim.rounds.roadsign.maps_file = 'roadsign_map.csv';

%% Process epochs
params.sim.epochs.difficulty.apply = true;
params.sim.epochs.difficulty.sequence_file = 'sequence_difficulty.csv';
params.sim.epochs.outcomes.apply = true;

%% Process events
params.sim.events.difficulty.apply = true;
params.sim.events.outcomes.apply = true;
params.sim.events.button_presses.apply = true;

% Includes traffic decision button presses?
params.sim.events.traffic_decision.apply = true;
params.sim.events.traffic_decision.write_to_file = true;

% Includes road sign changes?
params.sim.events.roadsign.apply = true;

% Includes saccades to road signs?
params.sim.events.roadsign_saccades.apply = true;



%% Plots
params.sim.plots.pass_color = [0.9 0.9 0.95];
params.sim.plots.baseline_color = [0.9 0.95 0.9];



