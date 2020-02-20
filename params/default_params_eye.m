%% This script defines default eye tracking parameters for the opends 
% processing pipeline

%% General stuff
params.eye.tracker_type = 'eyelink';
params.eye.sub_dir = 'eye';


%% Conversion
params.eye.convert.prefix = 'eye-';
params.eye.convert.exec_eyelink = 'edf2asc-mac';
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

% Plots
params.eye.blinks.plots.save = true;
params.eye.blinks.plots.offset = 10;


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

