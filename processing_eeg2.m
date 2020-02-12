%% Processes EEG data from OpenDS driving task
% 
% NOTE: If you want to remove bad ICA components, run "run_ica.m" before
% this script
%
% 1. Load raw or ICA-filtered data if available
% 3. Load and synchronize ET data - use to remove ocular artifacts from EEG
%    time series
% 4. Define trials from simulation events
% 5. Run Hilbert analysis
% 6. Run ERP analysis
% 7. Run time/frequency analysis
% 8. Summarize across subjects and perform statistical analyses
% 9. Plot results
%

%% Init - load parameters, bad channels, and subjects passing Q/C

clear;
load processing_eeg_params.mat
preproc = load('preproc_params_hd.mat');
proc = load('processing_params.mat');

addpath('../lib/fieldtrip-20190419/');
addpath('../lib/cPCOH/functions');
addpath('../lib/textprogressbar');

bad_channels = [];

show_plots = false;

if exist(params.eeg.bad_channel_file, 'file')
   opts = detectImportOptions(params.eeg.bad_channel_file);
   opts = setvartype(opts, 'char');
   bad_channels = readtable(params.eeg.bad_channel_file, opts);
   bad_channels.Properties.VariableNames = [{'Subject'},{'Channel'}];
end

load(proc.params.qc.file);
qc_score = cell2mat(qc_eeg(:,1));
idx = qc_score>1; %=proc.params.qc.cutoff;
subjects = qc_eeg(idx,2);

% Temp
params.eeg.erp.cfg.lpfreq = 15;
params.eeg.erp.cfg.hpfreq = 1;

params.eeg.erp.apply = false;
params.eeg.timefreq.apply = false;
params.eeg.hilbert.apply = false;
params.clobber = true;

params.eeg.erp.min_trials = 5;

%% Process all subjects

subjects = {'0739'};

for i = 1 : length(subjects)
    
    subject = subjects{i};
    
    outdir = sprintf('%s/%s/%s', preproc.params.root_dir, preproc.params.output_dir, subject);
    figdir = sprintf('%s/figures', outdir);
    
    flag_file = sprintf('%s/eeg_processing.done',outdir);
    
    if ~params.clobber && exist(flag_file, 'file')
        fprintf('\n== Subject %s already processed. Skipping. ==\n\n', subject);
        continue;
    end
    
    fprintf('\n== Starting processing for subject %s ==\n\n', subject);
    
    %% Load raw or ICA-filtered data if available

    % Loads EEG and Simulation Log events. If ICA-filtered data are available,
    % loads these instead. Time series for these events will be zeroed to the 
    % SimulatorStarted event.
    %
    
    data_file = sprintf('%s/processing_results_eeg_artfrej.mat',outdir);
    load( data_file );
    
    %% Load simulation events/epochs
    sim = get_simulation_events_eeg( outdir, preproc, proc, data );

    %% Define trials from simulation events
    
%     trials = get_trials_eeg( data, sim, params );


    %% Run Hilbert analysis

    if params.eeg.hilbert.apply
        
        fprintf('\nStarting Hilbert processing for %s...\n', subject);
        results = process_hilbert_eeg( data, params, preproc, outdir );
        save(sprintf('%s/processing_results_eeg_hilbert.mat',outdir), 'results', '-v7.3');
        fprintf('Done Hilbert processing for %s.\n', subject);
        
    end
    

    %% Run ERP analysis
    
    if params.eeg.erp.apply
        
        trials = get_trials_eeg( data, sim, params.eeg.erp );
        
        fprintf('\nStarting ERP processing for %s...\n', subject);
        results = process_erp_eeg( data, trials, params );
        save(sprintf('%s/processing_results_eeg_erp.mat',outdir), 'results', '-v7.3');
        fprintf('Done ERP processing for %s.\n', subject);
        
    end


    %% Run T/F analysis
    
    if params.eeg.timefreq.apply
        
        trials = get_trials_eeg( data, sim, params.eeg.timefreq );
        
        fprintf('\nStarting time/frequency processing for %s...\n', subject);
        results = process_timefreq_eeg( data, trials, params );
        save(sprintf('%s/%s.mat', outdir, params.eeg.timefreq.output_file), 'results', '-v7.3');
        fprintf('Done time/frequency processing for %s.\n', subject);
        
    end

    flag_file = sprintf('%s/eeg_processing.done',outdir);
    fid = fopen(flag_file, 'w+');
    fclose(fid);
    
    fprintf('\n== Finished processing for subject %s ==\n\n', subject);

end % Subject


