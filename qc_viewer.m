%% Manual QC by evaluating plots of the data

load preproc_params_hd.mat;

% subjects = csvimport(params.subject_list_file, 'noHeader', true, 'outputAsChar', true);
dat = fileread(params.subject_list_file);
subjects = strsplit(strtrim(dat),'\n');

subjects=[{'0739'}];
params.qc.view_only = true;

if ~params.qc.view_only
    if ~params.qc.clobber && exist(params.qc.file,'file')
        load(params.qc.file);
    else
        if exist(params.qc.file,'file'), delete(params.qc.file); end
        qc = cell(0,2);
    end
end

output_dir = sprintf('%s/%s', params.root_dir, params.output_dir);

%% For each subject
for i = 1 : length(subjects)
   
    close all;
    
    subject = subjects{i};
    results_dir = sprintf('%s/%s', output_dir, subject);
    results_file = sprintf('%s/results.mat',results_dir);
    
    if ~exist(results_file, 'file')
       fprintf('\nNo results found for subject %s. Skipping.\n', subject); 
    elseif ~params.qc.view_only && ~isempty(qc) && sum(ismember(qc(:,2),subject)) > 0
       fprintf('\nSubject %s already QC''d. Skipping.\n', subject);
    else
       fprintf('\nViewing subject %s.\n', subject); 
        
       load(results_file);
       
       h = plot_blinks ( results, data, params );
       resize_window(h,[1300 600],[200 800]);
       xlim([0 2]);
       
       h = plot_saccades2 ( results, data, params );
       resize_window(h,[1300 600],[200 800]);
       xlim([0 2]);
       
       h = plot_events2( results, params.plots.events, [{'Rounds'},{'LaneChanges'},{'Baseline'},{'SaccadeRate'},{'Overtakes'}] );
       resize_window(h,[1300 600],[200 800]);
       xlim([0 4]);
       
       if params.qc.view_only
           input('Press any key for next subject','s');
       else
           rating = input('Enter rating (neg to quit):');

           if rating<0
               fprintf('Terminated QC at subject %s', subject);
               break; 
           end

           % Insert into QC list
           qc = [qc(1:i-1,:);[{rating},{subject}];qc(i:end,:)];
           save(params.qc.file, 'qc');
       end
        
    end
    
    
    
end