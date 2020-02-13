%% Pre-process EEG data step 1 - Filtering and ICA

% This script implements a processing pipeline for EEG data. This data
% should already have been QC'd. Simulation/eye logs must already have been
% processed to run this pipeline ("preprocessing_eye.m").
%
% Note: you need to have loaded a variable named "params", with all the
% required fields, for this pipeline to run properly. Use the
% "set_preprocessing_params_eeg.m" script to set your parameters, then run
% this script.
%
% Steps:
% 1. Apply filters
%       a. Bandpass     [0.5-60Hz]
%       b. Notch        [50Hz]
% 2. Remove pre-defined bad channels
% 3. Run ICA on each participant
% 4. Save result as "eeg-{uid}-preproc1.mat"
%


%% 0. Read list of subjects and bad channel data

clear;

addpath('func');
addpath('plot');
addpath('../../lib/fieldtrip-20190419/');
addpath('../../lib/textprogressbar');

bad_channels = [];

load(params.general.qc.file);
qc_score = cell2mat(qc_eeg(:,1));
idx = qc_score>1;
subjects = qc_eeg(idx,2);


for i = 1 : length(subjects)
    
    subject = subjects{i};
    
    outdir = sprintf( '%s/%s', params.io.output_dir, subject );
    figdir = sprintf( '%s/figures', outdir );
    
    flag_file = sprintf( '%s/eeg_processing.done', outdir );
    
    if ~params.general.clobber && exist(flag_file, 'file')
        fprintf( '\n== Subject %s already processed. Skipping. ==\n\n', subject );
        continue;
    end
    
    if exist(flag_file, 'file')
       delete(flag_file); 
    end
    
    fprintf( '\n== Starting pre-processing step 1 for subject %s ==\n\n', subject );


    %% 1. Load data, synchronize with simulation/eye data
    
    % Load raw EEG data
    data = load_data_eeg( params, subject );
    if isempty(data)
        fprintf( '\n== Data acould not be loaded for subject %s. Skipping. ==\n\n', subject );
        continue;
    end
    
    % Synchronize with simlog/eye time series
    data = synchronize_data_eeg( params, data );

    fprintf( '\tDone loading data.\n' );
    
    
    %% 2. Remove bad channels
    data = remove_bad_channels_eeg( params, data );
    
    fprintf( '\tDone removing bad channels.\n' );
    
    
    %% 3. Apply bandpass+notch filters
    %       a. Bandpass     [default=0.5-60Hz]
    %       b. Notch        [default=50Hz]
    
    data = apply_bandpass_eeg( params, data );
    if isempty(data)
        fprintf( '\n== Could not apply filters for subject %s. Skipping. ==\n\n', subject );
        continue;
    end

    fprintf( '\tDone filtering.\n' );
    

    %% 4. Run ICA on each participant
    data = run_ica( params, data );
    if isempty(data)
        fprintf( '\n== Could not run ICA for subject %s. Skipping. ==\n\n', subject );
        continue;
    end
    
    fprintf( '\tDone ICA.\n' );
   

    %% 5. Save result
    output_file = sprintf( '%s/results_ica.mat', subj_dir );
    save( output_file, '-v7.3' );
    
    
    %% 6. Finalise
    flag_file = sprintf( '%s/ica.done', outdir );
    
    fprintf( '\n== Finished pre-processing step 1 for subject %s ==\n\n', subject );
    

end

