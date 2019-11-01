trialdefs = [{'all'},{'easy'},{'difficult'},{'positive'},{'negative'}];
events = [{'left_change'},{'right_change'}];

for i = 1 : length(subjects)
   
    subject = subjects{i};
    outdir = sprintf('%s/%s/%s', preproc.params.root_dir, preproc.params.output_dir, subject);
    results_file = sprintf('%s/processing_results_eeg_timefreq.mat',outdir);
    if ~ exist(results_file, 'file')
       warning('No results file found for subject %s!', subject);
       continue;
    end
    
    flag_file = sprintf('%s/timefreq_converted',outdir);
    
    if exist(flag_file, 'file')
       warning('Subject %s already converted. Skipping...', subject);
       continue;
    end
    
    T = load(results_file);
    
    results = [];
    
    for j = 1 : length(trialdefs)
        for e = 1 : length(events)
            result_j = T.results.eeg.timefreq.(events{e}).(trialdefs{j});
            new_result_j = [];
            new_result_j.trl = result_j.trl;
            new_result_j.trl_baseline = result_j.trl_baseline;
            new_result_j.channels = result_j.spectra{1}.channels;
            new_result_j.cfg = result_j.spectra{1}.cfg;
            new_result_j.spectra = {};

            for k = 1 : length(result_j.spectra)
                new_result_j.spectra(k) = {result_j.spectra{k}.stats.mean};
            end

            results.eeg.timefreq.(events{e}).(trialdefs{j}) = new_result_j;
        end
    end
    
    save(sprintf('%s/processing_results_eeg_timefreq.mat',outdir), 'results', '-v7.3');
    
    fclose(fopen(flag_file, 'w+'));

    fprintf('Done subject %s\n', subject);
    
end