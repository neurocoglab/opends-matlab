%% Pre-process EEG data step 3 - Artifact rejection
%
% This script implements a processing pipeline for EEG data. It identifies
% ocular and zvalue artifacts, removes these, and interpolates.
%
% It requires step 2 to have been completed ("preprocessing_2_eeg.m"), as
% well as the eye preprocessing and processing pipelines 
% ("preprocessing_eye.m" and "processing_eye.m").
%
% Note: you need to have loaded a variable named "params", with all the
% required fields, for this pipeline to run properly. Use the
% "default_params_X.m" scripts to set your parameters, then run
% this script.
%
% Steps:
% 1. Load results of step 2 (ICA component removal)
% 2. Synchronise time series
% 3. Identify saccade intervals as artifacts
% 4. Identify z-score outliers as artifacts
% 5. Remove artifacts and save result
%

%% 0. Load subject list

if ~exist('params','var')
    error('No params variable set. Aborting preprocessing.') 
end

addpath(params.general.fieldtrip_lib);

subjects = strsplit(strtrim(fileread(sprintf('%s/%s/%s', params.io.input_dir, ...
                                                         params.io.metadata_dir, ...
                                                         params.general.subjects_file))));

fprintf('\n\n==== START OF EEG PRE-PROCESSING STEP 3 ===\n\n');
                                             
fprintf('\nFound %d subjects.\n', length(subjects));


%% 1. For each subject, identify and remove artifacts

for i = 1 : length(subjects)
    
    subject = subjects{i};
    ok = true;
    
    clobber = params.general.clobber;

    outdir = sprintf( '%s/%s', params.io.output_dir, subject );
    figdir = sprintf( '%s/figures', outdir );
    flagdir = sprintf( '%s/flags', outdir );

    flag_done = sprintf('%s/eeg_preproc_2.done', flagdir);
    if ~exist(flag_done, 'file')
        fprintf('\n\tPreprocessing step 2 has not been completed for subject %s! Skipping.\n', subject);
        continue;
    end

    flag_file = sprintf( '%s/eeg_preproc_3.done', flagdir );

    if ~params.general.clobber && exist(flag_file, 'file')
        fprintf('\n\tPreprocessing step 3 already performed for subject %s! Skipping.\n', subject);
        continue;
    end

    if exist( flag_file, 'file' )
        delete( flag_file );
    end

    fprintf('\n\tProcessing subject %s...', subject);
    
    try

        if exist(flag_file, 'file')
            delete( flag_file );
        end

        input_file = sprintf('%s/results_preproc_eeg_2.mat', outdir);
        results_file = sprintf('%s/results_preproc_eeg_3.mat', outdir);
        load( input_file, 'data' );

        eye_file = sprintf('%s/results_preproc_eye.mat', outdir);
        T = load( eye_file );
        data.eye = T.data.eye;
        clear T;

        data = remove_artifacts_eeg( params, data );

        % Save result
        save(results_file, 'data', '-v7.3');

        % Plot?
        if params.eeg.artifacts.plots.save
            plot_artifacts_eeg( params, data, true );
        end

        fclose( fopen(flag_file,'w+') );

        close all;
        fprintf('done.\n');
    
    catch err
        if params.general.fail_on_error
            rethrow(err);
        end
        fprintf('\n== Errors encountered in pre-processing step 3 for subject %s: %s ==\n\n', subject, err.message);
        if params.general.show_error_stack
            fprintf(2, '%s\n', getReport(err, 'extended'));
        end
        
    end
    
end

fprintf('\n\n==== DONE EEG PRE-PROCESSING STEP 3 ===\n\n');

