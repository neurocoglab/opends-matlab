function [ summary ] = get_summary_eeg_erp( subjects, preproc )
%get_summary_eeg_erp Get summary of ERP processing for all subjects

added = false(length(subjects),1);
for i = 1 : length(subjects)
   
    subject = subjects{i};
    outdir = sprintf('%s/%s/%s', preproc.params.root_dir, preproc.params.output_dir, subject);
    flag_file = sprintf('%s/eeg_processing.done',outdir);
     if ~exist(flag_file, 'file')
        fprintf(' No ERP result for %s. Skipping.', subject);
     else
        added(i) = true;
     end
    
end

subs = subjects(added);

summary = [];

summary.erp.left_change.tlocked = cell(length(subs),1);
summary.erp.left_change.easy.tlocked = cell(length(subs),1);
summary.erp.left_change.difficult.tlocked = cell(length(subs),1);
summary.erp.left_change.positive.tlocked = cell(length(subs),1);
summary.erp.left_change.negative.tlocked = cell(length(subs),1);

summary.erp.right_change.tlocked = cell(length(subs),1);
summary.erp.right_change.easy.tlocked = cell(length(subs),1);
summary.erp.right_change.difficult.tlocked = cell(length(subs),1);
summary.erp.right_change.positive.tlocked = cell(length(subs),1);
summary.erp.right_change.negative.tlocked = cell(length(subs),1);

summary.erp.subjects = subs;

% fprintf(' Loading ERP results...');
textprogressbar(-1);
textprogressbar(sprintf(' Loading ERP results for %d subjects: ',length(subs)));

for i = 1 : length(subs)
   
    subject = subs{i};
    
    % Load ERP and add to summary
    outdir = sprintf('%s/%s/%s', preproc.params.root_dir, preproc.params.output_dir, subject);
    result_file = sprintf('%s/processing_results_eeg_erp.mat',outdir);
    T = load(result_file);
    
    summary.erp.left_change.tlocked(i) = {T.results.eeg.erp.left_change.timelock};
    summary.erp.left_change.easy.tlocked(i) = {T.results.eeg.erp.left_change.easy.timelock};
    summary.erp.left_change.difficult.tlocked(i) = {T.results.eeg.erp.left_change.difficult.timelock};
    summary.erp.left_change.positive.tlocked(i) = {T.results.eeg.erp.left_change.positive.timelock};
    summary.erp.left_change.negative.tlocked(i) = {T.results.eeg.erp.left_change.negative.timelock};

    summary.erp.right_change.tlocked(i) = {T.results.eeg.erp.right_change.timelock};
    summary.erp.right_change.easy.tlocked(i) = {T.results.eeg.erp.right_change.easy.timelock};
    summary.erp.right_change.difficult.tlocked(i) = {T.results.eeg.erp.right_change.difficult.timelock};
    summary.erp.right_change.positive.tlocked(i) = {T.results.eeg.erp.right_change.positive.timelock};
    summary.erp.right_change.negative.tlocked(i) = {T.results.eeg.erp.right_change.negative.timelock};
    
    textprogressbar(i/length(subs)*100);

end

textprogressbar('  done.');
% fprintf('done.\n');


end

