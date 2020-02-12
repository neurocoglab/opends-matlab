%% Set parameters and path

show_plots = false;

preproc = load('preproc_params_osx.mat');
load('processing_params_osx.mat');

% Apply QC filter
load(params.qc.file);
qc_score = cell2mat(qc(:,1));
idx = qc_score>=params.qc.cutoff;
subjects = qc(idx,2);
% subjects = [{'0133'}];

opt = detectImportOptions(params.subject_data_file);
opt.VariableOptions(1).Name = 'UID';

subject_data = readtable(params.subject_data_file, opt);

if ~exist(params.output_dir, 'dir')
   mkdir(params.output_dir); 
end

%% Define Screen ROIs

% This requires a simulation replay, producing screenshots locked to track
% positions.
% Screenshots will be labelled into:
%       


%% Load behavioural data, sync time series, label epochs


%% Load screenshots, label fixation regions



%% For baseline epochs, compute timelines, compute averages
%   Including:
%       Pupil diameter
%       Blink rate
%       Saccades (total)
%       Saccades to regions
%       Fixation time

clear summary;
summary.baseline.pupil = [];
summary.passing.pupil  = [];
summary.zscore.baseline.pupil = [];
summary.zscore.passing.pupil  = [];
summary.passing_diff.pupil  = [];
summary.zscore.passing_diff.pupil  = [];
summary.baseline.saccade_rate = [];
summary.passing.saccade_rate  = [];
summary.cycles.baseline.pupil = [];
summary.zscore.cycles.baseline.pupil = [];
summary.cycles.passing.pupil = [];
summary.zscore.cycles.passing.pupil = [];
summary.passing_outcome.positive.pupil = [];
summary.passing_outcome.negative.pupil = [];
summary.zscore.passing_outcome.positive.pupil = [];
summary.zscore.passing_outcome.negative.pupil = [];
summary.passing_diffs = {};
summary.passing_outcomes = {};



for i = 1 : length(subjects)
   
    subject = subjects{i};
    results_file = sprintf('%s/%s/results.mat',params.data_dir,subject);
    preprocess = load(results_file);
    
    results = []; % Load from file
    
    results.epochs = process_epochs( preprocess.results, preprocess.data, params, preproc );
    
    fn = sprintf('%s/%s/processing_results.mat', params.data_dir, subject);
    save(fn,'results');
    
    % Plots
    if params.show_plots
        plot_epoch_stats( results.epochs );

        % Save results and plots
    
    end
    
    % Aggregate summary stats
    summary.subjects = subjects;
    summary.diff_levels = results.epochs.diff_levels;
    summary.baseline.subjects.pupil(i) = {results.epochs.baseline.pupil};
    summary.passing.subjects.pupil(i) = {results.epochs.passing.pupil};
    summary.zscore.baseline.subjects.pupil(i) = {results.epochs.zscore.baseline.pupil};
    summary.zscore.passing.subjects.pupil(i) = {results.epochs.zscore.passing.pupil};
    
    summary.cycles.baseline.subjects.pupil(i) = {results.epochs.cycles.baseline.pupil};
    summary.cycles.passing.subjects.pupil(i) = {results.epochs.cycles.passing.pupil};
    summary.zscore.cycles.baseline.subjects.pupil(i) = {results.epochs.zscore.cycles.baseline.pupil};
    summary.zscore.cycles.passing.subjects.pupil(i) = {results.epochs.zscore.cycles.passing.pupil};
    
    summary.baseline.subjects.saccade_rate(i) = {results.epochs.baseline.saccade_rate};
    summary.passing.subjects.saccade_rate(i) = {results.epochs.passing.saccade_rate};
    
%     summary.baseline.pupil = [summary.baseline.pupil;results.epochs.baseline.pupil];
%     summary.passing.pupil = [summary.passing.pupil;results.epochs.passing.pupil];
%     summary.zscore.baseline.pupil = [summary.zscore.baseline.pupil;results.epochs.zscore.baseline.pupil];
%     summary.zscore.passing.pupil = [summary.zscore.passing.pupil;results.epochs.zscore.passing.pupil];
    
    N_cycles = length(preprocess.results.sim2track.cycle_times)+1;
    
    if isempty(summary.cycles.baseline.pupil)
        summary.cycles.baseline.pupil = results.epochs.cycles.baseline.pupil;
        summary.zscore.cycles.baseline.pupil = results.epochs.zscore.cycles.baseline.pupil;
        summary.cycles.passing.pupil = results.epochs.cycles.passing.pupil;
        summary.zscore.cycles.passing.pupil = results.epochs.zscore.cycles.passing.pupil;
    else
        N_cycles = min(N_cycles, length(summary.cycles.baseline.pupil));
        for j = 1 : N_cycles
            summary.cycles.baseline.pupil(j) = {[summary.cycles.baseline.pupil{j};results.epochs.cycles.baseline.pupil{j}]};
            summary.zscore.cycles.baseline.pupil(j) = {[summary.zscore.cycles.baseline.pupil{j};results.epochs.zscore.cycles.baseline.pupil{j}]};
            summary.cycles.passing.pupil(j) = {[summary.cycles.passing.pupil{j};results.epochs.cycles.passing.pupil{j}]};
            summary.zscore.cycles.passing.pupil(j) = {[summary.zscore.cycles.passing.pupil{j};results.epochs.zscore.cycles.passing.pupil{j}]};
        end
    end
    
    if isempty(summary.passing_diff.pupil)
        summary.passing_diff.pupil = results.epochs.passing_diff.pupil;
        summary.zscore.passing_diff.pupil = results.epochs.zscore.passing_diff.pupil;
    else
        for j = 1 : length(results.epochs.diff_levels)
            summary.passing_diff.pupil(j) = {[summary.passing_diff.pupil{j};results.epochs.passing_diff.pupil{j}]};
            summary.zscore.passing_diff.pupil(j) = {[summary.zscore.passing_diff.pupil{j};results.epochs.zscore.passing_diff.pupil{j}]};
        end
    end
    
    summary.idx_baseline(i) = {find(results.epochs.idx_baseline)};
    summary.idx_passing(i) = {find(results.epochs.idx_passing)};
    summary.passing_diffs(i) = {results.epochs.passing_diffs};
    summary.passing_outcomes(i) = {results.epochs.passing_outcomes};
    
    summary.passing_outcome.positive.subjects.pupil(i) = {results.epochs.passing_outcome.positive.pupil};
    summary.passing_outcome.negative.subjects.pupil(i) = {results.epochs.passing_outcome.negative.pupil};
    summary.zscore.passing_outcome.positive.subjects.pupil(i) = {results.epochs.zscore.passing_outcome.positive.pupil};
    summary.zscore.passing_outcome.negative.subjects.pupil(i) = {results.epochs.zscore.passing_outcome.negative.pupil};
    
    summary.baseline.saccade_rate = [summary.baseline.saccade_rate;results.epochs.baseline.saccade_rate];
    summary.passing.saccade_rate = [summary.passing.saccade_rate;results.epochs.passing.saccade_rate];
    
    fprintf('Finished epochs processing for subject %s\n', subject);
end

% Now do summary stats
summary.stats = analyze_epochs( summary, subject_data );

summary_file = sprintf('%s/%s', params.output_dir, params.epochs.summary_output);
save(summary_file , 'summary', '-v7.3' );

if params.epochs.save_plots
    plot_epoch_summary2( summary, params, true );
end

fprintf('Done epochs processing.\n');

%% Plot the above summary

% if show_plots
% 
%     load( summary_file );
%     h = plot_epoch_summary( summary, params, false, true );
% 
% end

%% For all events, compute pupil responses
%
%
%
%


%% Event-wise analysis
%
%   Choices following mistakes
%   Choices following collisions
%   Choices following correct choices
%   Early versus late
%   Experienced vs. novice
%   Good vs. poor performance
%   
%

ft_defaults;

summary = [];
summary.overtake.tlocked = cell(length(subjects),1);
summary.overtake.tlocked_bl = cell(length(subjects),1);
summary.overtake.tlocked_bl2 = cell(length(subjects),1);
summary.overtake.diffs = cell(length(subjects),1);
summary.overtake.outcomes = cell(length(subjects),1);

summary.left_change.tlocked = cell(length(subjects),1);
summary.left_change.tlocked_bl = cell(length(subjects),1);
summary.left_change.tlocked_bl2 = cell(length(subjects),1);
summary.left_change.diffs = cell(length(subjects),1);
summary.left_change.outcomes = cell(length(subjects),1);

summary.right_change.tlocked = cell(length(subjects),1);
summary.right_change.tlocked_bl = cell(length(subjects),1);
summary.right_change.tlocked_bl2 = cell(length(subjects),1);
summary.right_change.diffs = cell(length(subjects),1);
summary.right_change.outcomes = cell(length(subjects),1);

summary.subjects = subjects;

for i = 1 : length(subjects)
   
    subject = subjects{i};
    results_file = sprintf('%s/%s/results.mat',params.data_dir,subject);
    process_results_file = sprintf('%s/%s/processing_results.mat',params.data_dir,subject);
    preprocess = load(results_file);
    
    results = [];
    
    if exist(process_results_file,'file')
        load(process_results_file);
    end
    
    results.events = process_events( preprocess.results, preprocess.data, results, params, preproc );
    
%     x = results.events.overtake.tlocked;
    summary.overtake.tlocked(i) = {results.events.overtake.tlocked};
    summary.overtake.tlocked_bl(i) = {results.events.overtake.tlocked_bl};
    summary.overtake.tlocked_bl2(i) = {results.events.overtake.tlocked_bl2};
    summary.overtake.diffs(i) = {results.events.overtake.diffs};
    summary.overtake.outcomes(i) = {results.events.overtake.outcomes};

    summary.left_change.tlocked(i) = {results.events.left_change.tlocked};
    summary.left_change.tlocked_bl(i) = {results.events.left_change.tlocked_bl};
    summary.left_change.tlocked_bl2(i) = {results.events.left_change.tlocked_bl2};
    summary.left_change.diffs(i) = {results.events.left_change.diffs};
    summary.left_change.outcomes(i) = {results.events.left_change.outcomes};

    summary.right_change.tlocked(i) = {results.events.right_change.tlocked};
    summary.right_change.tlocked_bl(i) = {results.events.right_change.tlocked_bl};
    summary.right_change.tlocked_bl2(i) = {results.events.right_change.tlocked_bl2};
    summary.right_change.diffs(i) = {results.events.right_change.diffs};
    summary.right_change.outcomes(i) = {results.events.right_change.outcomes};
    
    if i == 1
        summary.overtake.t = results.events.overtake.t;
        summary.left_change.t = results.events.left_change.t;
        summary.right_change.t = results.events.right_change.t;
    end
    
    save(process_results_file,'results');
    
    fprintf('Finished events processing for subject %s\n', subject);
    
end

% Do summary stats
summary.stats = analyze_events( summary, params );

save( sprintf('%s/%s', params.output_dir, params.events.summary_output), 'summary', 'params', '-v7.3' );
fprintf('Done events processing.\n');


%% Plot the above summary

if show_plots

    load (sprintf('%s/%s', params.output_dir, params.events.summary_output));
    plot_event_stats2(summary, params);

end

%% 
