%% Print out notes for each subject with EEG

exp_dir = '/Volumes/AndrewElements/data/driving/exp';

load processing_eeg_params;
proc = load('processing_params.mat');

load(proc.params.qc.file);
qc_score = cell2mat(qc_eeg(:,1));
idx = qc_score>1; % =proc.params.qc.cutoff;
subjects = qc(idx,2);

output_file = 'eeg_notes.txt';
fidout = fopen(output_file, 'w');

qc_eeg = cell(0,2);

for s = 1 : length(subjects)
    subject = subjects{s};
    ok = true;
   
    %% Load data
    cfg = params.eeg.cfg;
    cfg.headerfile = sprintf('%s/%s/%s-eeg/%s.vhdr', params.eeg.data_dir, subject, subject, subject);
    if ~exist(cfg.headerfile, 'file')
       zip_file = sprintf('%s/%s/%s-eeg.zip', params.eeg.data_dir, subject, subject);
       if ~exist(zip_file, 'file')
           warning('Subject %s has no EEG data. Skipping...', subject);
           ok = false;
           qc_eeg = [qc_eeg;[{0},{subject}]];
       end
       
    end

    if ok
       
        qc_eeg = [qc_eeg;[{3},{subject}]];
        
        fprintf(fidout, '\n============ Subject: %s ============\n', subject);
        
        subj_dir = sprintf('%s/%s', exp_dir, subject);
        file = dir(sprintf('%s/*otes*', subj_dir));
        
        fid = fopen(sprintf('%s/%s', file.folder, file.name),'r');
        tline = fgets(fid);
        
        while tline > 0
            
            fprintf(fidout, tline);
            tline = fgets(fid);
        end
        
        fclose(fid);
        
        fprintf('Adding subject %s\n', subject);
        
    end
    
end

fclose(fidout);