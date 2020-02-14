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
% 1. Load data and synchronize with simulation/eye data
% 2. Remove pre-defined bad channels
% 3. Apply filters
%       a. Bandpass     [default=0.5-60Hz]
%       b. Notch        [default=50Hz]
% 4. Run ICA on each participant
% 5. Save result as "eeg-{uid}-preproc1.mat"
%

%% 0. Read list of subjects and bad channel data

if ~exist('params','var')
    error('No params variable set. Aborting preprocessing.') 
end

addpath(params.general.fieldtrip_lib);

subjects = strsplit(fileread(sprintf('%s/%s/%s', params.io.input_dir, ...
                                                 params.io.metadata_dir, ...
                                                 params.general.subjects_file)));

fprintf('\n\n==== START OF PROCESSING ===\n\n');
                                             
fprintf('\nFound %d subjects.\n', length(subjects));
                                 
for i = 1 : length(subjects)
    
    subject = subjects{i};
    
    outdir = sprintf( '%s/%s', params.io.output_dir, subject );
    figdir = sprintf( '%s/figures', outdir );
    
    % Check if output directory exists; if not, create it
    if ~exist(outdir, 'dir')
       mkdir(outdir); 
    end
    
    if ~exist(figdir, 'dir')
       mkdir(figdir); 
    end
    
    flag_file = sprintf( '%s/eeg_preproc_1.done', outdir );
    
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
    data = synchronize_data_eeg( params, data, subject );

    fprintf( '\tDone loading & synchronizing data.\n' );
    
    
    %% 2. Remove bad channels
    data = remove_bad_channels_eeg( params, data, subject );
    
    fprintf( '\tDone removing bad channels [%d found].\n', length(data.eeg.bad_channels) );
    
    
    %% 3. Apply bandpass+notch filters
    %       a. Bandpass     [default=0.5-60Hz]
    %       b. Notch        [default=50Hz]
    
    data = apply_filters_eeg( params, data );
    ftype = ''; conj = '';
    if params.eeg.bandpass.apply
        ftype = 'low/high-pass ';
        conj = 'and ';
    end
    if params.eeg.notch.apply
        ftype = [ftype conj 'notch '];
    end

    fprintf( '\tDone applying %sfiltering.\n', ftype );
    

    %% 4. Run ICA on each participant
    
    fprintf( '\tRunning ICA for %s... (may take a while)\n', subject );
    data = run_ica( params, data );
    if isempty(data)
        fprintf( '\n== Could not run ICA for subject %s. Skipping. ==\n\n', subject );
        continue;
    end
    
    fprintf( '\tDone ICA.\n' );
   

    %% 5. Save result
    output_file = sprintf( '%s/results_preproc_eeg_1.mat', outdir );
    save( output_file, '-v7.3' );
    
    
    %% 6. Finalise
    fclose( fopen( flag_file, 'w' ) );
    
    fprintf( '\n== Finished pre-processing step 1 for subject %s ==\n\n', subject );
    

end

fprintf('\n\n==== END OF PROCESSING ===\n\n');

