%% Pre-process EEG data step 2 - Interactive ICA component rejection
%
% This script allows the user to interactively select and reject individual
% ICA components. It requires step 1 to have been completed
% ("preprocessing_1_eeg.m").
%
% Note: you need to have loaded a variable named "params", with all the
% required fields, for this pipeline to run properly. Use the
% "set_preprocessing_params_eeg.m" script to set your parameters, then run
% this script.
%
% Steps:
% 1. Load ICA result
% 2. Show ICA results
% 3. Prompt for components to reject
% 4. Remove components and save result as "eeg-{uid}_preproc2.mat"
%
%

if ~exist('params','var')
    error('No params variable set. Aborting preprocessing.') 
end

addpath(params.general.fieldtrip_lib);

subjects = strsplit(strtrim(fileread(sprintf('%s/%s/%s', params.io.input_dir, ...
                                                         params.io.metadata_dir, ...
                                                         params.general.subjects_file))));

fprintf('\n\n==== START OF PROCESSING ===\n\n');
                                             
fprintf('\nFound %d subjects.\n', length(subjects));
       
%% Loop through subjects and prompt for input

for i = 1 : length(subjects)
    
    try 
        
        subject = subjects{i};
        ok = true;

        clobber = params.general.clobber;

        outdir = sprintf( '%s/%s', params.io.output_dir, subject );
        figdir = sprintf( '%s/figures', outdir );
        flagdir = sprintf( '%s/flags', outdir );

        flag_done = sprintf( '%s/eeg_preproc_1.done', flagdir );

        if ~exist(flag_done, 'file')
            fprintf('\n\tSubject %s has no ICA results! Skipping.\n', subject);
            continue;
        end

        flag_file = sprintf( '%s/eeg_preproc_2.done', flagdir );

        if exist(flag_file, 'file') && ~clobber
            fprintf('\n\tICA rejection already performed for subject %s! Skipping.\n', subject);
            continue;
        end

        fprintf('\n\tProcessing subject %s.\n', subject);

        if exist(flag_file, 'file')
            delete( flag_file );
        end

        input_file = sprintf('%s/results_preproc_eeg_1.mat', outdir);
        results_file = sprintf('%s/results_preproc_eeg_2.mat', outdir);
        load( input_file, 'data' );

        cfg = params.eeg.ica.plots.cfg;
        h = figure;
        h.Color = 'w';
        evalc('ft_topoplotIC(cfg, data.eeg.ica)');
        resize_window(h, params.eeg.ica.plots.window_size_topo);

        cfg.viewmode = 'component';
        evalc('ft_databrowser(cfg, data.eeg.ica)');
        pause(1);
        resize_window(gcf, params.eeg.ica.plots.window_size_browser);

        % Prompt for bad components
        channels = input(sprintf('\tSubject %s: Enter ICA components to remove (space delimited): ',subject), 's');
        to_rem = cellfun(@str2num, strsplit(channels, ' '));

        if ~isempty(to_rem)
           fprintf('\n\tRemoving %d components...', length(to_rem));
           cfg = [];
           cfg.component = to_rem;
           [~,data.eeg.ft] = evalc('ft_rejectcomponent(cfg, data.eeg.ica, data.eeg.ft)');

        end

        data.eeg.ica = [];
        data.eeg.ica.removed = to_rem;

        % Save result
        save(results_file, 'data', '-v7.3');
        if params.eeg.ica.plots.save
            saveas(gcf,sprintf('%s/ica_browser_eeg.fig',figdir));
            saveas(h,sprintf('%s/ica_topoplot_eeg.png',figdir));
        end

        fclose( fopen(flag_file,'w+') );

        close all;
        fprintf('done.\n');
    
    catch err
        fprintf('\n== Errors encountered in pre-processing step 2 for subject %s: %s ==\n\n', subject, err.message);
    end

end