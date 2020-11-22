%% This script defines default eye tracking parameters for the opends 
% processing pipeline

%% General stuff
params.eye.sub_dir = 'eye';


%% Conversion
params.eye.convert.format = 'eyelink';
params.eye.convert.prefix = 'eye-';
params.eye.convert.exec_eyelink = 'edf2asc.exe';
params.eye.convert.exec_smi = 'smi2csv';
params.eye.convert.columns = 'smi_columns.txt';

%% Eye tracking
params.eye.Fs = 500;

%% Gap detection
params.eye.gaps.apply = true;
params.eye.gaps.criterion = 'diam';
params.eye.gaps.gapmin = 10;
params.eye.gaps.gapthres = 5;
params.eye.gaps.gapmerge = 400;
params.eye.gaps.plots.color = [0.9 0.8 0.8];
params.eye.gaps.buffer = 0;

%% Blink detection
params.eye.blinks.apply = true;
params.eye.blinks.criterion = 'diam';
params.eye.blinks.smooth = true;
params.eye.blinks.smooth_width = 1000;
params.eye.blinks.gapbuffer = 50;
params.eye.blinks.window = 10000;
params.eye.blinks.absthres = 10;
params.eye.blinks.interval = [70,120];
params.eye.blinks.thres = [-0.5 0.1];
params.eye.blinks.maxblink = 100;
params.eye.blinks.from_gaps.width_lims = [];

% Plots
params.eye.blinks.plots.save = true;
params.eye.blinks.plots.offset = 10;
params.eye.blinks.plots.color = [0.8 0.9 0.8];


%% Saccade detection
params.eye.saccades.apply = true;
params.eye.saccades.xpos_variable = 'pos_x';
params.eye.saccades.ypos_variable = 'pos_y';
params.eye.saccades.monitor_dims = [530,300];
params.eye.saccades.monitor_res = [1920,1080];
params.eye.saccades.distance = 700;
params.eye.saccades.nitr_smooth = 10;
params.eye.saccades.velocity_thres = 40;
params.eye.saccades.min_width = 0;
params.eye.saccades.max_width = 1000;
params.eye.saccades.sr_window = 12000;
params.eye.saccades.peaks_file = 'saccade_peak.mat';

% Plots
params.eye.saccades.plots.save = true;
params.eye.saccades.plots.show_lines = true;
params.eye.saccades.plots.show_rate = true;


%% Luminance correction
params.eye.luminance.apply = true;
params.eye.luminance.sub_dir = 'luminance';
params.eye.luminance.downsample = 5;
params.eye.luminance.smooth = 0;
params.eye.luminance.use_offset = 0;
params.eye.luminance.offsets = -2:0.1:2;
params.eye.luminance.screen_pct = 50;
params.eye.luminance.outlier_lim = 6;

% Plots

%% Processing - epochs
params.eye.epochs.apply = true;
params.eye.epochs.sub_dir = 'epochs';
params.eye.epochs.overtake_window = [-10 10];
params.eye.epochs.plots.save = true;
params.eye.epochs.plots.zlims = [-4 4];
params.eye.epochs.plots.scatter = false;
params.eye.epochs.plots.show_webplots = true;



%% Processing - events
params.eye.events.zscore = true;
params.eye.events.smooth = 500;

params.eye.events.difficulty.apply = true;
params.eye.events.outcomes.apply = true;

params.eye.events.plots.save = true;
params.eye.events.plots.ylims = [-1.5 1];
params.eye.events.plots.show_webplots = true;

params.eye.events.overtake.zscore = true;
params.eye.events.overtake.prepost = [5000 5000];
% params.eye.events.overtake.baseline = [2000 500];
params.eye.events.overtake.plots.ylims = [-1 1.5];

params.eye.events.left_change.zscore = true;
params.eye.events.left_change.prepost = [7500 5000];
% params.eye.events.left_change.baseline = [10000 8000];
params.eye.events.left_change.plots.ylims = [-1 1.5];

params.eye.events.right_change.zscore = true;
params.eye.events.right_change.prepost = [5000 7500];
% params.eye.events.right_change.baseline = [2000 500];
params.eye.events.right_change.plots.ylims = [-1 1.5];

params.eye.events.saccades.zscore = true;
params.eye.events.saccades.vmin = 100;
params.eye.events.saccades.prepost = [1500 1500];
% params.eye.events.saccades.baseline = [2000 500];
params.eye.events.saccades.plots.ylims = [-0.5 0.5];

params.eye.events.blinks.zscore = true;
params.eye.events.blinks.prepost = [1500 1500];
% params.eye.events.blinks.baseline = [2000 500];
params.eye.events.blinks.plots.ylims = [-0.5 0.5];

params.eye.events.random.N = 10;

params.eye.events.alpha = 0.05;
params.eye.events.min_trials = 5;


