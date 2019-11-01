%% Unzip

for s = 1 : length(subjects)
    subject = subjects{s};
    ok = true;

    %% Load data
    cfg = params.eeg.cfg;
    cfg.headerfile = sprintf('%s/%s/%s-eeg/%s.vhdr', params.eeg.data_dir, subject, subject, subject);
    if ~exist(cfg.headerfile, 'file')
       zip_file = sprintf('%s/%s/%s-eeg.zip', params.eeg.data_dir, subject, subject);
       unzip_dir = sprintf('%s/%s/%s-eeg', params.eeg.data_dir, subject, subject);
       fprintf('Zip file: %s\n', zip_file);
       if ~exist(zip_file, 'file')
           warning('Subject %s has no EEG data. Skipping...', subject);
           ok = false;
       else
           unzip(zip_file, unzip_dir);
           ok = exist(cfg.headerfile, 'file');
       end
    end
    
end

%% Clean up

for s = 1 : length(subjects)
    subject = subjects{s};
    ok = true;
   
    unzip_dir = sprintf('%s/%s/%s-eeg', params.eeg.data_dir, subject, subject);
    zip_file = sprintf('%s/%s/%s-eeg.zip', params.eeg.data_dir, subject, subject);
    
    if exist(zip_file, 'file')
        if exist(unzip_dir, 'dir')
            rmdir(unzip_dir, 's');
            fprintf('Removed unzip data for %s\n', subject);
        end
    end
    
    
    
end


%% Add flag files

for s = 1 : length(subjects)
    subject = subjects{s};
    outdir = sprintf('%s/%s/%s', preproc.params.root_dir, preproc.params.output_dir, subject);
    flag_file = sprintf('%s/eeg-erp.done', outdir);
    result_file = sprintf('%s/processing_results_eeg.mat',outdir);
    if exist(result_file, 'file')
        fid = fopen(flag_file, 'w');
        fclose(fid);
        fprintf('Done subject %s\n', subject);
    else
        fprintf('Skipped subject %s\n', subject);
    end
    
end
