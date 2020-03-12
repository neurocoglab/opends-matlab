%% Processes EEG data from OpenDS driving task
% 
% NOTE: Run this script only after running "preprocessing_eeg.m"
%       A "params" variable must have already been set in the environment;
%       see the "default_params_X.m" scripts in opends-project
% 
%
% 1. Load preprocessed data
% 2. Run Hilbert analysis
% 3. Run ERP analysis
% 4. Run time/frequency analysis
% 5. Summarize across subjects and perform statistical analyses
% 6. Plot results
%

%% Load subject list

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

clobber = params.general.clobber;

fprintf('\n\n==== START OF EEG PROCESSING ===\n\n');
                                             
fprintf('\nFound %d subjects.\n', length(subjects));


%% Run Hilbert processing

if params.eeg.hilbert.apply
    
    flag_file = sprintf('%s/hilbert_eeg.done', results_flagdir);

    if exist(flag_file, 'file')
       delete( flag_file ); 
    end
    
    summary = [];

    for i = 1 : length( subjects )
        
        try

            subject = subjects{i};
        
            fprintf('\nStarting Hilbert processing for %s...\n', subject);
            
            outdir = sprintf( '%s/%s', params.io.output_dir, subject );
            figdir = sprintf( '%s/figures', outdir );
            flagdir = sprintf( '%s/flags', outdir );
            
            results_file = sprintf('%s/results_hilbert_eeg.mat', outdir);

            flag_file_i = sprintf('%s/eeg_hilbert.done', flagdir);
            if exist(flag_file_i, 'file') && ~clobber
               fprintf('\tHilbert analysis already performed for subject %s.\n', subject);
               continue;
            end

            flag_done = sprintf('%s/eeg_preproc_3.done', flagdir);
            if ~exist(flag_done, 'file')
               fprintf('\tEEG preprocessing step 3 not run for subject %s! Skipping...\n', subject);
               continue;
            end
            
            if exist(flag_file_i, 'file'), delete(flag_file_i); end
            input_file = sprintf('%s/results_preproc_eeg_3.mat', outdir);
            load( input_file );
            
            results = process_hilbert_eeg( params, data, [] );
            
            save(results_file, 'results', '-v7.3');
            
            if params.eeg.hilbert.plots.save
                load( sprintf('%s/results_preproc_eye.mat', outdir) );
                plot_hilbert_eeg( params, data, results, true );
            end
            
            fclose( fopen( flag_file_i, 'w+' ) );
            
            clear data results;
            
            fprintf('done.\n');
        
        catch err
            warning on;
            warning('\nError encountered while processing %s:\n%s\n', subject, err.message);

        end

    end
    
    fprintf('\nDone Hilbert processing.\n');
    
end

%% Run statistical analyses on Hilbert output

if params.eeg.hilbert.apply
    
    summary = analyze_hilbert_eeg( params );

    summary_file = sprintf('%s/summary_hilbert_eeg.mat', results_dir);
    save(summary_file , 'summary', '-v7.3' );

    if params.eeg.hilbert.plots.save
        plot_hilbert_summary_eeg( params, summary, true );
    end

    fclose( fopen( flag_file, 'w+' ) );

    fprintf('\nDone Hilbert analysis.\n');
        
end
    

%% Run ERP analysis

if params.eeg.erp.apply

    flag_file = sprintf('%s/erp_eeg.done', results_flagdir);

    if exist(flag_file, 'file')
       delete( flag_file ); 
    end
    
    summary = [];

    for i = 1 : length( subjects )
        
%         try
        
            fprintf('\nStarting ERP processing for %s...\n', subject);
            
            outdir = sprintf( '%s/%s', params.io.output_dir, subject );
            figdir = sprintf( '%s/figures', outdir );
            flagdir = sprintf( '%s/flags', outdir );
            
            results_file = sprintf('%s/results_eeg.mat', outdir);

            flag_file_i = sprintf('%s/eeg_erp.done', flagdir);
            if exist(flag_file_i, 'file') && ~clobber
               fprintf('\tERP analysis already performed for subject %s.\n', subject);
               load( results_file );
               summary = update_erp_summary( params, results, summary );
               continue;
            end

            flag_done = sprintf('%s/eeg_preproc_3.done', flagdir);
            if ~exist(flag_done, 'file')
               fprintf('\tEEG preprocessing step 3 not run for subject %s! Skipping...\n', subject);
               continue;
            end
            
            if exist(flag_file_i, 'file'), delete(flag_file_i); end
            
            input_file = sprintf('%s/results_preproc_eeg_3.mat', outdir);
            load( input_file );
            
            results = [];
            
            load( results_file );
            
            [results, summary] = process_erp_eeg( params, data, results, summary );
            
            save(results_file, 'results', '-v7.3');
            
            if params.eeg.erp.plots.save
                plot_erp_eeg( params, results, true );
            end
            
            clear data results;
            
            fclose( fopen( flag_file_i, 'w+' ) );
            
            fprintf('done.\n');
        
%         catch err
%             warning on;
%             warning('\nError encountered while processing %s:\n%s\n', subject, err.message);
% 
%         end

    end
    
    % Run statistical analyses
    summary = analyze_erp_eeg( params, summary );

    summary_file = sprintf('%s/summary_erp_eeg.mat', results_dir);
    save(summary_file , 'summary', '-v7.3' );

    if params.eeg.erp.plots.save
        plot_erp_summary_eeg( params, summary, true );
    end

    fclose( fopen( flag_file, 'w+' ) );

    fprintf('\nDone ERP processing & analysis.\n');

end


%% Run T/F analysis

if params.eeg.timefreq.apply

    flag_file = sprintf('%s/timefreq_eeg.done', results_flagdir);

    if exist(flag_file, 'file')
       delete( flag_file ); 
    end
    
    summary = [];

    for i = 1 : length( subjects )
        
%         try
        
            fprintf('\nStarting time/frequency processing for %s...\n', subject);
            
            outdir = sprintf( '%s/%s', params.io.output_dir, subject );
            figdir = sprintf( '%s/figures', outdir );
            flagdir = sprintf( '%s/flags', outdir );

            flag_file_i = sprintf('%s/eeg_timefreq.done', flagdir);
            if exist(flag_file_i, 'file') && ~clobber
               fprintf('\tTime/freq analysis already performed for subject %s.\n', subject);
               continue;
            end

            flag_done = sprintf('%s/eeg_preproc_3.done', flagdir);
            if ~exist(flag_done, 'file')
               fprintf('\tTime/frequency preprocessing step 3 not run for subject %s! Skipping...\n', subject);
               continue;
            end
            
            if exist(flag_file_i, 'file'), delete(flag_file_i); end
            
            input_file = sprintf('%s/results_preproc_eeg_3.mat', outdir);
            load( input_file );
            
            results = [];
            results_file = sprintf('%s/results_eeg.mat', outdir);
            load( results_file );
            
            [results] = process_timefreq_eeg( params, data, results );
            
            save(sprintf('%s/results_eeg.mat', outdir), 'results', '-v7.3');
            
            if params.eeg.timefreq.plots.save
                plot_timefreq_eeg( params, results, true );
            end
            
            clear data results;
            
            fclose( fopen( flag_file_i, 'w+' ) );
            
            fprintf('done.\n');
        
%         catch err
%             warning on;
%             warning('\nError encountered while processing %s:\n%s\n', subject, err.message);
% 
%         end

    end
    
    fprintf('\nDone time/frequency processing.\n');
    
end

    
%% Run statistical analyses

if params.eeg.timefreq.apply

    fprintf('\nStarting time/frequency analysis\n');
    
    summary = analyze_timefreq_eeg( params, summary );

    summary_file = sprintf('%s/summary_timefreq_eeg.mat', results_dir);
    save(summary_file , 'summary', '-v7.3' );

    if params.eeg.timefreq.plots.save
        plot_timefreq_summary_eeg( params, summary, true );
    end

    fclose( fopen( flag_file, 'w+' ) );

    fprintf('\nDone time/frequency analysis.\n');
    
end



fprintf('\n\n==== DONE EEG PROCESSING ===\n\n');

