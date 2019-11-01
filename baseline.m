%% Baseline analysis

% Author 2017 Andrew Reid
% Requires preprocess.m to already have been executed

%% 0. Initialize



%% 1. Read baseline distance intervals

load('baseline_intervals.mat');
temp = load('preproc_params.mat');
params.preproc = temp.params;
temp = load('baseline_params.mat');
params.baseline = temp.params.baseline;
clear temp;

%% 2. For each subject (average)

%   2.1. Read preprocessed pupil time series
pupil = load(sprintf('%s/blink.mat', params.preproc.output_dir));

%   2.2. Read simulator event log from CSV files


%   2.3. Match distances to time points (interpolate where necessary)
%   2.4. Extract baseline pupil time series
%   2.5. Average across each interval
%   2.6. Plot time series and averages

%% 3. For each subject (fixations)
%   3.1. Compute saccades and fixation events
%   3.2. Obtain fixation-onset-locked epochs
%   3.3. Obtain average ER pupil responses for baseline
%   3.4. Plot average ERPRs

