%% This script defines default simulation log parameters for the opends 
% processing pipeline

%% General
params.sim.sub_dir = 'driving';

%% Convert log files
params.sim.convert.prefix = 'simlog-';
params.sim.convert.filter = 'serial-filter';
params.sim.convert.exec = 'simlog2csv';
