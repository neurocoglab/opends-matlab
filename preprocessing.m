%% Eye tracking data preprocessing steps

%% 1. Load parameters and subject list

clear;

load preproc_params_osx.mat;

% Temp DEBUG
params.clobber=true;

% subjects = csvimport(params.subject_list_file, 'noHeader', true, 'outputAsChar', true);
dat = fileread(params.subject_list_file);
subjects = strsplit(strtrim(dat),'\n');

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
        
        outdir = sprintf('%s/%s/%s', params.root_dir, params.output_dir, subject);
        if ~exist(outdir,'dir')
           mkdir(outdir);
        elseif params.clear_existing
           fprintf('Clearing existing data...\n');
%            rmdir(outdir,'s');
%            mkdir(outdir);
        end

        results_file = sprintf('%s/results.mat',outdir);

        flag_file = sprintf('%s/eye_logs.done', outdir);

        if ~exist(flag_file,'file') || clobber

           clobber = true;  % If this is redone, the rest must also be redone
            
           fprintf('Converting eye tracking log for %s...', subject);

           if exist(flag_file, 'file'), delete(flag_file); end 

           rdir = params.root_dir;
           if ~isempty(params.data_dir)
              rdir = sprintf('%s/%s', rdir, params.data_dir);
           end
           
           output_dir = sprintf('%s/%s/%s', rdir, ...
                                               params.convert.input_dir, ...
                                               subject);

           log_file = sprintf('%s/%s_samples.txt', output_dir, ...
                                                   subject);
           
           log_files = [];
           parts = [];
                                               
           if ~exist(log_file, 'file')
                [a,fn,ext] = fileparts(log_file);
                log_file = sprintf('%s/%s_part1%s', a, fn, ext);
                k = 2;
                while exist(log_file, 'file')
                    log_files = [log_files {log_file}];
                    log_file = sprintf('%s/%s_part%d%s', a, fn, k, ext);
                    k = k + 1;
                end
           else
                log_files = {log_file};
           end
           
           if isempty(log_files)
              
              log_file = sprintf('%s/%s_samples.txt', output_dir, ...
                                                      subject);
                                               
              zip_file = sprintf('%s/%s%s.zip', output_dir, ...
                                                params.convert.prefix, ...
                                                subject);
      
              try
                  if exist(zip_file, 'file')
                     % Extract samples log from zip file
                     unzip(zip_file, output_dir);
                     
                     % Is is in the current directory?
                     lfile = dir(sprintf('%s/*Samples.txt', output_dir));
                     
                     if ~isempty(lfile)
                         lfile = lfile(1);
                         movefile(sprintf('%s/%s', output_dir, lfile.name), ...
                                  sprintf('%s/tmp', output_dir));
                         delete(sprintf('%s/*.txt', output_dir));
                         movefile(sprintf('%s/tmp', output_dir), ...
                                  log_file);
                         log_files = {log_file};
                     else
                         % Otherwise, find subdirectory containing subject
                         % ID (if this doesn't exist we have a problem)
                         subdirs = dir(output_dir);
                         subdirs = subdirs(3:end);
                         subdirs = subdirs([subdirs.isdir]);
                         subdir = [];
                         for s = 1 : length(subdirs)
                            if strfind(subdirs(s).name, subject)
                               subdir = subdirs(s).name;
                            end
                         end
                         if ~isempty(subdir)
                             lfiles = dir(sprintf('%s/%s/*Samples.txt', output_dir, subdir));
                        
                             for f = 1 : length(lfiles)
                                 lfile = lfiles(f);
                                 fname = lfile.name;
                                 idx = strfind(fname,'part');
                                 if idx > 0
                                     k = fname(idx+4);
                                     lfile2 = sprintf('%s/%s/%s', output_dir, subdir, fname);
                                     out_file = log_file;
                                     if length(lfiles) > 1
                                         [a,fn,ext] = fileparts(out_file);
                                         out_file = sprintf('%s/%s_part%s%s', a, fn, k, ext);
                                     end
                                     movefile(lfile2, out_file);
                                     log_files = [log_files {out_file}];
                                     parts = [parts {k}];
%                                      fprintf('Created %s\n', out_file);
                                 else
                                    warning('No part in sample file name: %s', fname)
                                 end
                             end
                             rmdir(sprintf('%s/%s/%s', output_dir, subject),'s');
                         end
                     end
                  end
              catch err
                  % This will get handled by the next statement 
                  a=0;
              end
           end

           if isempty(log_files)
               fprintf('\nLog or zip file not found for %s; skipping subject.\n', subject);
               ok = false;
           else
               prefix = subject;
               for j = 1 : length(log_files)
                   log_file = log_files{j};
                   if length(log_files) > 1
                       idx = strfind(log_file,'part');
                       %fprintf('%s\n',log_file);
                       if idx > 0
                           k = log_file(idx+4);
                       end
                       prefix = sprintf('%s_part%s', subject, k);
                   end
                   
                   a = strfind(params.convert.exec,'/');
                   a = a(end);
                   cdir = params.convert.exec(1:a);
                   cfile = params.convert.exec(a+1:end);
                   cmd = sprintf('cd "%s"; ./%s "%s" -o "%s" -columns "%s" -prefix %s -clobber', ...
                                    cdir, ...
                                    cfile, ...
                                    log_file, ...
                                    outdir, ...
                                    params.convert.columns, ...
                                    prefix);
%                     fprintf('%s\n',cmd);
                    [status,result] = system(cmd);
                    if status ~= 0
                        fprintf('\nError converting log for %s.\n%s', subject, result);
                        ok = false;
                    else
                        fid = fopen(flag_file,'w');
                        fclose(fid);
                    end
               end
               if ok
                  fprintf('Done.\n');
               else
                  fprintf('Done with errors.\n');
               end
           end
        else
            fprintf('Eye tracking log for "%s" already converted; skipping.\n', subject);
        end

        
        
        %% Load converted time series
        if ok
            
            input_files = [];
            input_file = sprintf('%s/%s/%s/%s_samples.csv', params.root_dir, params.output_dir, subject, subject);
            
            if exist(input_file, 'file')
               input_files = {input_file};
            else
               input_file = sprintf('%s/%s/%s/%s_part1_samples.csv', params.root_dir, params.output_dir, subject, subject);
               k = 2;
               while exist(input_file, 'file')
                   input_files = [input_files {input_file}];
                   input_file = sprintf('%s/%s/%s/%s_part%d_samples.csv', params.root_dir, params.output_dir, subject, subject, k);
                   k = k + 1;
               end
            end
            
            flag_file = sprintf('%s/%s/%s/time_series.done', params.root_dir, params.output_dir, subject);

            if exist(flag_file,'file') && ~clobber
%                 load(results_file, 'data');
                fprintf('Times series for %s already converted; skipping.\n', subject);
            else

                clobber = true;  % If this is redone, the rest must also be redone
                
                fprintf('Converting time series for %s...', subject');

                if exist(flag_file,'file'), delete(flag_file); end
                
                t_start = NaN;
                Fs = NaN; ts = NaN;
                data = []; data_i = [];
                data.subject = subject;
                data_i.subject = subject;

                for j = 1 : length(input_files)
                   
                    input_file = input_files{j};
                    T = import_samples(input_file);
                    
                    data_i.eye.t = T{1}(2:end); % cell2mat(T{1}(2:end));
                    data_i.eye.pos_left_x = T{3}(2:end); % cell2mat(T(2:end,3));
                    data_i.eye.pos_left_y = T{4}(2:end); % cell2mat(T(2:end,4));
                    data_i.eye.diam_left = T{5}(2:end); % cell2mat(T(2:end,5));

                    clear T;

                    % Convert double to single
                    data_i.eye.pos_left_x = single(data_i.eye.pos_left_x);
                    data_i.eye.pos_left_y = single(data_i.eye.pos_left_y);
                    data_i.eye.diam_left = single(data_i.eye.diam_left);

                    % Convert time to ms relative to start
                    data_i.eye.t = data_i.eye.t / 1000;
                    if j == 1
                        data.eye.t_start = data_i.eye.t(1);
                        ts = double(data_i.eye.t(2)-data_i.eye.t(1));
                        Fs = params.Fs;
                        ts = 1/Fs*1000;
                        data_i.eye.t_start = data.eye.t_start;
                    else
                        data_i.eye.t_start = data.eye.t_start + data_i.eye.t(1);
                    end
                    
                    data_i.eye.t = data_i.eye.t - data.eye.t_start;
                    data_i.eye.t = single(data_i.eye.t);
                    data_i.eye.Fs = Fs;

                    if j == 1
                       data = data_i;
                    else
                       % Append to existing data
%                        data_i.eye.tgap(:,1) = data_i.eye.tgap(:,1) + length(data.eye.t); % Offset indices
                       data.eye.t = [data.eye.t;data_i.eye.t];
%                        data.eye.tgap = [data.eye.tgap;data_i.eye.tgap];
                       data.eye.pos_left_x = [data.eye.pos_left_x;data_i.eye.pos_left_x];
                       data.eye.pos_left_y = [data.eye.pos_left_y;data_i.eye.pos_left_y];
                       data.eye.diam_left = [data.eye.diam_left;data_i.eye.diam_left];
                    end
                    
                end
               
                % Find gaps in sampling
                data.eye = remove_gaps2( data.eye, params.blink );
                
                save(results_file,'data');
                
                % Only set flag if there is a single log; otherwise
                % there is still work to do...
                if length(input_files) == 1
                    fid = fopen(flag_file,'w');
                    fclose(fid);
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

            flag_file = sprintf('%s/%s/%s/sim_logs.done', params.root_dir, params.output_dir, subject);
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
                if status ~= 0 || ~isempty(strfind(result,'Exception'))
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
            results = process_saccades2( results, data.eye, params.saccades );
            plot_saccades2 ( results, data, params, true );

            % Save results
            save(results_file, 'data', 'results');
            fid = fopen(flag_file,'w');
            fclose(fid);
            fprintf('Done.\n');
        end

        if params.show_plots
            plot_saccades2 ( results, data, params );
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
            
            load(results_file);
            
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