function [ summary ] = get_summary_eeg_timefreq( subjects, preproc, params )
%get_summary_eeg_timefreq Get summary of time/frequency processing for all subjects

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

summary.timefreq.left_change.all.tlocked = cell(length(subs),1);
summary.timefreq.left_change.easy.tlocked = cell(length(subs),1);
summary.timefreq.left_change.difficult.tlocked = cell(length(subs),1);
summary.timefreq.left_change.positive.tlocked = cell(length(subs),1);
summary.timefreq.left_change.negative.tlocked = cell(length(subs),1);

summary.timefreq.right_change.all.tlocked = cell(length(subs),1);
summary.timefreq.right_change.easy.tlocked = cell(length(subs),1);
summary.timefreq.right_change.difficult.tlocked = cell(length(subs),1);
summary.timefreq.right_change.positive.tlocked = cell(length(subs),1);
summary.timefreq.right_change.negative.tlocked = cell(length(subs),1);

summary.timefreq.subjects = subs;

% fprintf(' Loading ERP results...');
textprogressbar(-1);
textprogressbar(sprintf(' Loading time/frequency results for %d subjects: ',length(subs)));

% has_inf = false(length(subs),8);

for i = 1 : length(subs)
   
    subject = subs{i};
    
    % Load ERP and add to summary
    outdir = sprintf('%s/%s/%s', preproc.params.root_dir, preproc.params.output_dir, subject);
    result_file = sprintf('%s/%s.mat',outdir,params.eeg.timefreq.output_file);
    T = load(result_file);

    summary.timefreq.left_change.all.tlocked(i) = {T.results.eeg.timefreq.left_change.all.timelock};
    summary.timefreq.left_change.easy.tlocked(i) = {T.results.eeg.timefreq.left_change.easy.timelock};
    summary.timefreq.left_change.difficult.tlocked(i) = {T.results.eeg.timefreq.left_change.difficult.timelock};
    summary.timefreq.left_change.positive.tlocked(i) = {T.results.eeg.timefreq.left_change.positive.timelock};
    summary.timefreq.left_change.negative.tlocked(i) = {T.results.eeg.timefreq.left_change.negative.timelock};

    summary.timefreq.right_change.all.tlocked(i) = {T.results.eeg.timefreq.right_change.all.timelock};
    summary.timefreq.right_change.easy.tlocked(i) = {T.results.eeg.timefreq.right_change.easy.timelock};
    summary.timefreq.right_change.difficult.tlocked(i) = {T.results.eeg.timefreq.right_change.difficult.timelock};
    summary.timefreq.right_change.positive.tlocked(i) = {T.results.eeg.timefreq.right_change.positive.timelock};
    summary.timefreq.right_change.negative.tlocked(i) = {T.results.eeg.timefreq.right_change.negative.timelock};
    
    textprogressbar(i/length(subs)*100);

end

textprogressbar('  done.');
% fprintf('done.\n');


end

