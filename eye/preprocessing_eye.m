%% Eye tracking / simulation log preprocessing steps
%
% This script performs pre-processing on eye and simulation log data. A
% "params" variable must already have been set when this script is run. Use
% "default_params_X.m" to set default parameters.
%
% 1. Convert eye tracking log from native/export format to CSV format
% 2. Load converted time series and save as a Matlab struct
% 3. Convert the simulation log (XML format) to CSV files
% 4. Detect and interpolate over eye blinks
% 5. Identify saccade and fixation intervals
% 6. Regress out estimated relative screen luminance signal (requires
%    gaze-related luminance to already have been estimated)
% 7. Identify round and repeat intervals from simulation log and save as a
%    Matlab struct
%

%% Load subject list

if ~exist('params','var')
    error('No params variable set. Aborting preprocessing.') 
end

addpath(params.general.fieldtrip_lib);

subjects = strsplit(strtrim(fileread(sprintf('%s/%s/%s', params.io.input_dir, ...
                                                         params.io.metadata_dir, ...
                                                         params.general.subjects_file))));

                                                     
fprintf('\n\n==== START OF EYE/SIM PRE-PROCESSING ===\n\n');

fprintf('\nFound %d subjects.\n', length(subjects));

%% For each subject
for i = 1 : length(subjects)
    
    subject = subjects{i};
    ok = true;
    
    clobber = params.general.clobber;

    fprintf('\n-- Processing subject %s --\n\n', subject);
    
    outdir = sprintf( '%s/%s', params.io.output_dir, subject );
    figdir = sprintf( '%s/figures', outdir );
    flagdir = sprintf( '%s/flags', outdir );
    
    try
        
        %% 1. Convert log to CSV
       
        % Check if output directory exists; if not, create it
        if ~exist(outdir, 'dir')
           mkdir(outdir); 
        end

        if ~exist(figdir, 'dir')
           mkdir(figdir); 
        end
        
        if ~exist(flagdir, 'dir')
           mkdir(flagdir); 
        end

        results_file = sprintf('%s/results_preproc_eye.mat',outdir);
        flag = 'eye_logs_converted.done';
        flag_file = sprintf('%s/%s', flagdir, flag);

        if ~exist(flag_file,'file') || clobber

           clobber = true;  % If this is redone, the rest must also be redone
            
           fprintf('\tConverting %s log for %s...', params.eye.convert.format, subject);
           delete_flags( flag, flagdir );

           warning on;
           
           switch params.eye.convert.format
            
               case 'smi'
                   ok = convert_eyelog_smi( subject, params );
                   
               case 'eyelink'
                   ok = convert_eyelog_eyelink( subject, params );

               case 'tobii'
                   ok = convert_eyelog_tobii( subject, params );
                   
               otherwise
                   % Fail outright
                   error('\tEye tracker type "%s" is unsupported!', params.eye.convert.format)

           end
          
           if ok 
              fclose( fopen(flag_file,'w') );
              fprintf('done.\n'); 
           else
              fprintf('\n');
              warning('Errors encountered.'); 
           end
           
        else
            fprintf('\tEye tracking log for "%s" already converted; skipping.\n', subject);
        end

        
        
        %% 2. Load converted time series
        if ok
            
            flag = 'eye_logs_imported.done';
            flag_file = sprintf('%s/%s', flagdir, flag);

            if exist(flag_file,'file') && ~clobber
                fprintf('\tTimes series for %s already converted; skipping.\n', subject);
            else

                clobber = true;  % If this is redone, the rest must also be redone
                
                fprintf('\tLoading converted time series for %s...', subject');

                delete_flags( flag, flagdir );
                
                data = read_log_eye( params, subject );
                
                if ~isempty(data)
                    save(results_file, 'data');
                    fclose( fopen(flag_file,'w') );
                    clear data;
                    fprintf(' done.\n');
                else
                    ok = false;
                    warning('Could not read converted time series. Skipping subject %s.', subject);
                end
             
            end
            
        end
       
        
        
    %% 3. Convert simulation event log to CSV files & synchronize time series
        if ok
      
            flag = 'sim_logs_converted.done';
            flag_file = sprintf('%s/%s', flagdir, flag);
            
            if exist(flag_file, 'file') && ~clobber
                fprintf('\tSim log for %s already converted; skipping.\n', subject);
            else
                
                clobber = true;  % If this is redone, the rest must also be redone
                fprintf('\tConverting sim log for %s...', subject);

                delete_flags( flag, flagdir );
                
                [ok] = convert_log_sim( params, subject );
                
                if ok
                    load( results_file );
                    data = synchronize_eye_sim( params, data );
                    save(results_file, 'data');
                    fclose( fopen(flag_file,'w') );
                    fprintf(' done.\n');
                else
                    ok = false;
                    warning('Could not convert simulation log. Skipping subject %s.', subject);
                end
                
            end

        end
        

    %% 4. Eye blink removal

        if ok

            flag = 'eye_blinks.done';
            flag_file = sprintf('%s/%s', flagdir, flag);

            if exist(flag_file, 'file') && ~clobber
                fprintf('\tEye blinks for %s already processed; skipping.\n', subject);
            else
                delete_flags( flag, flagdir );
                
                clobber = true;  % If this is redone, the rest must also be redone

                load(results_file);

                fprintf('\tProcessing eyeblinks for %s...', subject);
                data = process_blinks_eye( data, params );
                
                % Save results
                save(results_file, 'data');
                if params.eye.blinks.plots.save
                    plot_blinks_eye ( data, params, true );
                end
                
                fclose( fopen(flag_file,'w') );
                fprintf(' done.\n');
            end

            if params.general.show_plots
                plot_blinks_eye ( data, params, false );
            end

            clear data;

        end                           

           
    %% 5. Identify saccade and fixation intervals

    if ok

        flag = 'eye_saccades.done';
        flag_file = sprintf('%s/%s', flagdir, flag);

        if exist(flag_file, 'file') && ~clobber
            fprintf('\tSaccades for %s already processed; skipping.\n', subject);
        else

            delete_flags( flag, flagdir );
            
            clobber = true;  % If this is redone, the rest must also be redone

            load(results_file);

            fprintf('\tProcessing saccades for %s...', subject);
            data = process_saccades_eye( data, params );
            
            if params.eye.saccades.plots.save
                plot_saccades_eye ( data, params, true );
            end

            % Save results
            save(results_file, 'data');
            fclose( fopen(flag_file,'w') );
            fprintf(' done.\n');
        end

        if params.general.show_plots
            plot_saccades_eye ( data, params, false );
        end

        clear data;
        
    end                           

    %% 6. Regress out estimated relative screen luminance signal
    
    if ok
        
        flag = 'eye_luminance.done';
        flag_file = sprintf('%s/%s', flagdir, flag);
        
        if exist(flag_file, 'file') && ~clobber
            fprintf('\tLuminance correction for %s already done; skipping.\n', subject);
        else

            fprintf('\tProcessing luminance for %s...', subject);
            
            load(results_file);
            
            warning off;
            data = process_luminance_eye( data, params );
            warning on;

            if ~isfield(data.eye, 'luminance')
                warning(' No luminance data found... skipping!');
            else
                clobber = true;
                delete_flags( flag, flagdir );

                if data.eye.luminance.deficient
                    warning(' Regression was rank-deficient... skipping!');
                    plot_luminance_eye ( data, params, true );
                else
                    plot_luminance_eye ( data, params, true );
                end

                % Save results
                save(results_file, 'data');
                fclose( fopen(flag_file,'w') );
                fprintf(' done.\n');
            
                if params.general.show_plots
                    plot_luminance_eye ( data, params, false );
                end
                
                clear data;
                
            end

        end
        
    end
    
    
    %% 7. Identify round and repeat intervals

    if ok
        
        flag = 'sim_rounds.done';
        flag_file = sprintf('%s/%s', flagdir, flag);

        if exist(flag_file, 'file') && ~clobber
            fprintf('\tRounds + baseline for %s already processed; skipping.\n', subject);
        else

            delete_flags( flag, flagdir );
            
            clobber = true;  % If this is redone, the rest must also be redone

            load(results_file);

            fprintf('\tProcessing rounds for %s...', subject);
            data = process_rounds_sim( data, params );

            % Save results and delete previous
            save(results_file, 'data');
            
            if params.sim.rounds.plots.save
                plot_events_sim( params, data, [{'Rounds'},{'LaneChanges'},{'Baseline'},{'SaccadeRate'},{'Overtakes'}], true );
            end

            fclose( fopen(flag_file,'w') );
            fprintf(' done.\n');

            % Plot rounds 
            if params.general.show_plots
                h = plot_events_sim( params, data, [{'Rounds'},{'LaneChanges'},{'Baseline'},{'SaccadeRate'},{'Overtakes'}] );
            end
            
            clear data;

        end

    end
    
    %% 8. Process fixation to areas of interest (AOIs)

    if ok && params.sim.events.fixations.apply
        
        flag = 'sim_fixation_events.done';
        flag_file = sprintf('%s/%s', flagdir, flag);

        if exist(flag_file, 'file') && ~clobber
            fprintf('\tFixation events for %s already processed; skipping.\n', subject);
        else

            delete_flags( flag, flagdir );
            
            clobber = true;  % If this is redone, the rest must also be redone

            load(results_file);

            fprintf('\tProcessing fixation events for %s...', subject);
            data = process_fixation_events( data, params );

            % Save results and delete previous
            save(results_file, 'data');
            
            if params.sim.rounds.plots.save
                plot_events_sim( params, data, [{'Rounds'},{'LaneChanges'},{'Baseline'},{'SaccadeRate'},{'Overtakes'}], true );
            end

            fclose( fopen(flag_file,'w') );
            fprintf(' done.\n');

            % Plot rounds 
            if params.general.show_plots
                h = plot_events_sim( params, data, [{'Rounds'},{'LaneChanges'},{'Baseline'},{'SaccadeRate'},{'Overtakes'}] );
            end
            
            clear data;

        end

    end

    catch err
        if params.general.show_error_stack
            disp(getReport(err, 'extended'));
           % fprintf(2, '%s\n', getReport(err, 'extended'));
        else
            warning on;
            warning('\nError encountered while processing %s:\n%s\n', subject, err.message);
        end
        ok = false;
    end
    
    if ok
        fprintf('\n-- Done subject %s --\n\n', subject);
    else
        fprintf('\n-- Failed for subject %s --\n\n', subject);
        if params.general.fail_on_error
            break;
        end
    end

    close all;

end

fprintf('\n\n==== DONE EYE/SIM PRE-PROCESSING ===\n\n');


%%
% Delete this flags and all flags after it
function delete_flags ( start_flag, flagdir )

    flag_files = [{'eye_logs_converted.done'}, ...
                  {'eye_logs_imported.done'}, ...
                  {'sim_logs_converted.done'}, ...
                  {'eye_blinks.done'}, ...
                  {'eye_saccades.done'}, ...
                  {'eye_luminance.done'}, ...
                  {'sim_rounds.done'}, ...
                  {'sim_fixation_events.done'}];
              
    idx = find(strcmp(flag_files, start_flag));
    if ~isempty(idx)
       for i = idx : length(flag_files)
           flag_file = sprintf('%s/%s', flagdir, flag_files{i});
           if exist(flag_file, 'file'), delete(flag_file); end
       end
    end
    


end