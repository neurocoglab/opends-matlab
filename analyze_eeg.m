%% Analyze EEG results
%
% Requires that processing_eeg.m has already been run

%% Init

clear;

addpath('../lib/fieldtrip-20190419');
addpath('../lib_areid/');
addpath('../lib/textprogressbar');
addpath('../lib/bwperimtrace');
addpath('../lib/plotly/plotly');

load processing_eeg_params.mat
preproc = load('preproc_params_hd.mat');
proc = load('processing_params.mat');

show_plots = true;
if ~exist(params.eeg.output_dir, 'dir')
   mkdir(params.eeg.output_dir); 
end
figures_dir = [params.eeg.output_dir '/figures'];
if ~exist(figures_dir, 'dir')
   mkdir(figures_dir); 
end

load(proc.params.qc.file);
qc_score = cell2mat(qc_eeg(:,1));
idx = qc_score>1; %=proc.params.qc.cutoff;
subjects = qc_eeg(idx,2);


%% Epoch-based
%
% Analyse Hilbert envelope over epochs of interest across subjects
%
%




%% Event-based (ERP & Time/frequency)
%
% Analyze ERP & time/frequency across subjects
%
%

show_plots = false;

summary = get_summary_eeg_erp( subjects, preproc );
plot_erp_subjects(summary, 'Fz', [-0.5 1], [-0.2 0.3]);
stats = analyze_eeg_erp( summary, params );
plot_eeg_erp_stats( stats, params, figures_dir, show_plots );
save(sprintf('%s/stats_eeg_erp.mat', params.eeg.output_dir), 'summary', 'stats');
clear summary stats;

summary = get_summary_eeg_timefreq( subjects, preproc, params );
stats = analyze_eeg_timefreq( summary, params, sprintf('%s/tmp_stats_timefreq', params.eeg.output_dir) );
plot_eeg_timefreq_stats( stats, params, figures_dir, show_plots );
save(sprintf('%s/stats_eeg_timefreq.mat', params.eeg.output_dir), 'stats','-v7.3');
% save(sprintf('%s/stats_eeg_timefreq.mat', params.eeg.output_dir), 'summary', 'stats','-v7.3');
clear summary stats;



