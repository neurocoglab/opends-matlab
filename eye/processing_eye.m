%% Eye tracking processing steps
%
% This script performs processing and statistical analysis on eye tracking
% data. This requires that all pre-processing steps already be completed.
%
%

if ~exist('params','var')
    error('No params variable set. Aborting preprocessing.') 
end

addpath(params.general.fieldtrip_lib);
ft_defaults;

subjects = strsplit(strtrim(fileread(sprintf('%s/%s/%s', params.io.input_dir, ...
                                                         params.io.metadata_dir, ...
                                                         params.general.subjects_file))));

metadata_file = sprintf('%s/%s/%s', params.io.input_dir, ...
                                    params.io.metadata_dir, ...
                                    params.general.subject_metadata_file);

opt = detectImportOptions(metadata_file);
opt.VariableOptions(1).Name = 'UID';
subject_data = readtable(metadata_file, opt);

results_dir = params.io.results_dir;
if ~exist(results_dir, 'dir')
   mkdir( results_dir ); 
end

results_figdir = sprintf('%s/figures', results_dir);
if ~exist(results_figdir, 'dir')
   mkdir( results_figdir ); 
end

results_flagdir = sprintf('%s/flags', results_dir);
if ~exist(results_flagdir, 'dir')
   mkdir( results_flagdir ); 
end
                                                     
fprintf('\n\n==== START OF PROCESSING ===\n\n');
                                             
fprintf('\nFound %d subjects.\n', length(subjects));

clobber = params.general.clobber;

ok = true;

%% Process epochs
%
%  Compute averages
%  Compare classes (baseline, passing, difficulty, outcome)
%  Compile summary
%  Run statistical analyses
%  Plot results
%

if params.eye.epochs.apply

    flag_file = sprintf('%s/epochs_eye.done', results_flagdir);

    if exist(flag_file, 'file')
       delete( flag_file ); 
    end

    summary = [];

    fprintf('\nProcessing epochs...\n');

    for i = 1 : length(subjects)

        try

            subject = subjects{i};

            outdir = sprintf( '%s/%s', params.io.output_dir, subject );
            figdir = sprintf( '%s/figures', outdir );
            flagdir = sprintf( '%s/flags', outdir );

            flag_file_i = sprintf('%s/eye_epochs.done', flagdir);
            if exist(flag_file_i, 'file') && ~clobber
               fprintf('\tEpochs already processed for subject %s! Skipping...\n', subject);
               continue;
            end

            flag_done = sprintf('%s/sim_rounds.done', flagdir);
            if ~exist(flag_done, 'file')
               fprintf('\tRounds not processed for subject %s! Skipping...\n', subject);
               continue;
            end

            if exist(flag_file_i, 'file'), delete(flag_file_i); end

            input_file = sprintf('%s/results_preproc_eye.mat',outdir);
            load(input_file);

            [results, summary] = process_epochs_eye( params, data, summary );

            results_file = sprintf('%s/results_eye.mat', outdir);
            save(results_file, 'results');

            % Plots
            if params.eye.epochs.plots.save
                plot_epoch_subject_eye( results, params, true );
            end

            fclose( fopen( flag_file_i, 'w+' ) );

            fprintf('\tFinished subject %s\n', subject);

        catch err
            warning on;
            warning('\nError encountered while processing %s:\n%s\n', subject, err.message);

        end

    end

    % Run statistical analyses
    summary = analyze_epochs_eye( summary, subject_data );

    summary_file = sprintf('%s/summary_epochs_eye.mat', results_dir);
    save(summary_file , 'summary', '-v7.3' );

    if params.eye.epochs.plots.save
        plot_epoch_summary_eye( params, summary, true );
    end

    fclose( fopen( flag_file, 'w+' ) );

    fprintf('\nDone epochs processing.\n');

end


%% Process events
%
%  Compute time-locked averages
%  Compare classes (baseline, passing, difficulty, outcome)
%  Compile summary
%  Run statistical analyses
%  Plot results
%

if params.eye.events.apply

    flag_file = sprintf('%s/events_eye.done', results_flagdir);


    if exist(flag_file, 'file')
       delete( flag_file ); 
    end

    summary = [];

    fprintf('\nProcessing events...\n');

    for i = 1 : length(subjects)

        try
            subject = subjects{i};

            outdir = sprintf( '%s/%s', params.io.output_dir, subject );
            figdir = sprintf( '%s/figures', outdir );
            flagdir = sprintf( '%s/flags', outdir );

            flag_file_i = sprintf('%s/eye_events.done', flagdir);
            if exist(flag_file_i, 'file') && ~clobber
               fprintf('\tEvents already processed for subject %s! Skipping...\n', subject);
               continue;
            end

            flag_done = sprintf('%s/sim_rounds.done', flagdir);
            if ~exist(flag_done, 'file')
               fprintf('\tRounds not processed for subject %s! Skipping...\n', subject);
               continue;
            end

            if exist(flag_file_i, 'file'), delete(flag_file_i); end

            input_file = sprintf('%s/results_preproc_eye.mat',outdir);
            load(input_file);

            results_file = sprintf('%s/results_eye.mat', outdir);
            if exist(results_file, 'file')
                load(results_file);
            else
                results = [];
            end

            [results, summary] = process_events_eye( params, data, results, summary );

            save(results_file, 'results');

            fclose( fopen( flag_file_i, 'w+' ) );

            fprintf('\tFinished subject %s\n', subject);

        catch err
            if params.general.show_error_stack
                fprintf(2, '%s\n', getReport(err, 'extended'));
                if params.general.fail_on_error
                    ok = false;
                    break;
                end
            else
                warning on;
                warning('\nError encountered while processing %s:\n%s\n', subject, err.message);
            end
        end

    end


    %% Run statistical analyses
    if ok
        summary = analyze_events_eye( params, summary );

        summary_file = sprintf('%s/summary_events_eye.mat', results_dir);
        save(summary_file , 'summary', '-v7.3' );

        if params.eye.events.plots.save
            plot_event_summary_eye( params, summary, true );
        end

        fclose( fopen( flag_file, 'w+' ) );

        fprintf('\nDone events processing.\n');
    end

    
end


fprintf('\n\n==== DONE PROCESSING ===\n\n');

