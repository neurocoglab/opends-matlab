%% Eye tracking data preprocessing steps

%% 1. Load parameters and subject list

if ~exist('params','var')
    error('No params variable set. Aborting preprocessing.') 
end

addpath('../driving2019');
addpath('../lib/fieldtrip-20190419');

% load preproc_params.mat;

% subjects = csvimport(params.subject_list_file, 'noHeader', true, 'outputAsChar', true);
opts = detectImportOptions(params.subject_list_file);
subject_table = readtable(params.subject_list_file, opts);
subjects = subject_table.Subject;

fprintf('\n\n==== START OF PROCESSING ===\n\n');

%% For each subject
for i = 1 : length(subjects)
    
    clear results;

    subject = subjects{i};
    ok = true;
    
    clobber = params.clobber;

    fprintf('\n-- Processing subject %s --\n\n', subject);
    
    try
        
        %% Convert log to CSV
        output_dir = sprintf('%s/%/%ss', params.root_dir, params.data_dir, params.output_dir);
        if ~exist(output_dir,'dir')
           mkdir(output_dir); 
        end
        
        outdir = sprintf('%s/%s/%s/%s', params.root_dir, params.data_dir, params.output_dir, subject);
        if ~exist(outdir,'dir')
           mkdir(outdir);
        elseif params.clear_existing
           fprintf('Clearing existing data...\n');
        end

        results_file = sprintf('%s/results.mat',outdir);
        flag_file = sprintf('%s/eye_logs.done', outdir);

        if ~exist(flag_file,'file') || clobber

           clobber = true;  % If this is redone, the rest must also be redone
            
           fprintf(' Converting %s log for %s...', params.tracker_type, subject);

           if exist(flag_file, 'file'), delete(flag_file); end 

           switch params.tracker_type
            
               case 'smi'
                   ok = convert_eyelog_smi( subject, params );
                   
               case 'eyelink'
                   ok = convert_eyelog_eyelink( subject, params );
                   
               otherwise
                   % Fail outright
                   error(' Eye tracker type "%s" is invalid!', params.tracker_type)
                   
           end
           
           if ok 
              fid = fopen(flag_file,'w'); fclose(fid);
              fprintf(' Done.\n'); 
           else
              fprintf(' Errors encountered.\n'); 
           end
           
        else
            fprintf(' Eye tracking log for "%s" already converted; skipping.\n', subject);
        end

        
        
        %% Load converted time series
        if ok
            
            flag_file = sprintf('%s/%s/%s/%s/time_series.done', params.root_dir, params.data_dir, params.output_dir, subject);

            if exist(flag_file,'file') && ~clobber
                fprintf('Times series for %s already converted; skipping.\n', subject);
            else

                clobber = true;  % If this is redone, the rest must also be redone
                
                fprintf('Converting time series for %s...', subject');

                if exist(flag_file,'file'), delete(flag_file); end
                
                data = read_eye_data( subject, params );
                
                if ~isempty(data)
                    save(results_file,'data');
                    fid = fopen(flag_file,'w');
                    fclose(fid);
                    clear data;
                end
                
                fprintf('Done.\n');

            end
        end
       
        
        
    %% Convert simulation event log to CSV files
        if ok
            
            data_dir = params.root_dir;
            if ~isempty(params.data_dir)
               data_dir = sprintf('%s/%s', data_dir, params.data_dir); 
            end

            flag_file = sprintf('%s/%s/%s/sim_logs.done', data_dir, params.output_dir, subject);
            if exist(flag_file, 'file') && ~clobber

                fprintf('Sim log for %s already converted; skipping.\n', subject);

            else
                
                clobber = true;  % If this is redone, the rest must also be redone

                fprintf('Converting sim log for %s...', subject);

                if exist(flag_file, 'file'), delete(flag_file); end

                if ~isempty(params.convert_sim.log_dir)
                   log_dir = sprintf('%s/%s', params.convert_sim.log_dir, subject);
                else
                   log_dir = subject;
                end
                log_file = sprintf('%s/%s/%s%s.log', data_dir, log_dir, params.convert_sim.prefix, subject);

                a = strfind(params.convert_sim.exec,'/');
                a = a(end);
                cdir = params.convert_sim.exec(1:a);
                cfile = params.convert_sim.exec(a+1:end);

                cmd = sprintf('cd "%s"; ./%s "%s" assets/opends/%s.xml "%s" -addserial', cdir, ...
                                                            cfile, ...
                                                            log_file, ...
                                                            params.convert_sim.filter, ...
                                                            outdir);

                [status,result] = system(cmd);
                if status ~= 0 || contains(result,'Exception')
                    fprintf('\nError converting simlog for %s.\n%s\n', subject, result);
                    ok = false;
                else
                    fid = fopen(flag_file,'w');
                    fclose(fid);
                    fprintf('Done.\n');
                end                               

            end

        end
        

    %% 2. Eye blink removal

        if ok

            flag_file = sprintf('%s/%s/%s/blinks.done', params.root_dir, params.output_dir, subject);

            if exist(flag_file, 'file') && ~clobber
                fprintf('Eye blinks for %s already processed; skipping.\n', subject);
%                 load(results_file);
            else
                if exist(flag_file, 'file'), delete(flag_file); end
                
                clobber = true;  % If this is redone, the rest must also be redone

                load(results_file);

                fprintf('Processing eyeblinks for %s...', subject);
                results = process_eyeblinks( data.eye, params.blink );
                
                results.subject = subject;
                results.params = params;

                % Save results
                save(results_file, 'data', 'results');
                plot_blinks ( results, data, params, true );

                fid = fopen(flag_file,'w');
                fclose(fid);
                fprintf('Done.\n');
            end

            if params.show_plots
                plot_blinks ( results, data, params );
            end

            %clear results data;

        end                           

    %% 4. Correction for foreshortening effect
    
    
    
    
        
    %% 5. Identify saccade and fixation intervals

    if ok

    %     output_file = sprintf('%s/saccades.mat',params.output_dir);
        flag_file = sprintf('%s/%s/%s/saccades.done', params.root_dir, params.output_dir, subject);

        if exist(flag_file, 'file') && ~clobber
            fprintf('Saccades for %s already processed; skipping.\n', subject);
%             load(results_file);
        else

            if exist(flag_file, 'file'), delete(flag_file); end
            
            clobber = true;  % If this is redone, the rest must also be redone

            load(results_file);

            fprintf('Processing saccades for %s...', subject);
            results = process_saccades( results, data.eye, params.saccades );
            plot_saccades ( results, data, params, true );

            % Save results
            save(results_file, 'data', 'results');
            fid = fopen(flag_file,'w');
            fclose(fid);
            fprintf('Done.\n');
        end

        if params.show_plots
            plot_saccades ( results, data, params );
        end

    end                           

    %% 6. Regress out estimated relative screen luminance signal
    
    if ok
        
        flag_file = sprintf('%s/%s/%s/luminance.done', params.root_dir, params.output_dir, subject);
        
        if exist(flag_file, 'file') && ~clobber
            fprintf('Luminance correction for %s already done; skipping.\n', subject);
%             load(results_file);
        else
            
            fprintf('Processing luminance for %s...', subject);
            
            warning off;
            
            lum_file = sprintf('%s/%s/%s/luminance.csv', params.root_dir, params.luminance.dir, data.subject);
            if ~exist(lum_file, 'file')
                warning on;
                warning(' No luminance values found... skipping!');
            else
            
                results = process_luminance( results, data, params );
                warning on;
                
                if results.luminance.deficient
                    warning(' Regression was rank-deficient... skipping!');
                    plot_luminance ( results, data, params, true );
                else
                    plot_luminance ( results, data, params, true );
                end
                
                % Save results
                save(results_file, 'data', 'results');
                fid = fopen(flag_file,'w');
                    fclose(fid);
                fprintf('Done.\n');
            end
            
        end
    
        if params.show_plots
            plot_luminance ( results, data, params );
        end
        
    end
    
    
    %% 7. Identify round and repeat intervals

    % TODO: use SimulationStart and SimulationEnd events in new version of
    % ar.OpenDS

    if ok

        %output_file = sprintf('%s/%s/sim2track.mat', params.data_dir, subject);
        flag_file = sprintf('%s/%s/%s/sim2track.done', params.root_dir, params.output_dir, subject);

        if exist(flag_file, 'file') && ~clobber

            fprintf('Rounds + baseline for %s already processed; skipping.\n', subject);
%             load(results_file);

        else

            if exist(flag_file, 'file'), delete(flag_file); end
            
            clobber = true;  % If this is redone, the rest must also be redone

            load(results_file);

            fprintf('Processing rounds for %s...', subject);
            [data, results] = process_rounds( data, results, params, subject );

            % Save results and delete previous
            save(results_file, 'data', 'results');
            plot_events( results, params, data, [{'Rounds'},{'LaneChanges'},{'Baseline'},{'SaccadeRate'},{'Overtakes'}], true );

            fid = fopen(flag_file,'w');
            fclose(fid);
            fprintf('Done.\n');

            % Plot rounds 
            if params.show_plots
                h = plot_events( results, params, data, [{'Rounds'},{'LaneChanges'},{'Baseline'},{'SaccadeRate'},{'Overtakes'}] );
            end

        end

    end

    catch err
        warning on;
        warning('\nError encountered while processing %s:\n%s\n', subject, err.message);
        ok = false;
    end
    
    if ok
        fprintf('\n-- Done subject %s --\n\n', subject);
    else
        fprintf('\n-- Failed for subject %s --\n\n', subject);
    end

    close all;

end

fprintf('\n\n==== END OF PROCESSING ===\n\n');

