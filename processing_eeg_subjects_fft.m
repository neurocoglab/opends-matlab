
addpath C:\Users\marij\Documents\MATLAB\Toolboxes\fieldtrip-20180918
ft_defaults

%% params

% preproc params
load processing_eeg_params.mat
preproc = load('preproc_params_hd.mat');
proc = load('processing_params.mat');

% get rid of this
dirm = 'C:\Users\marij\Documents\Projects_Driving\';
preproc.params.root_dir = [dirm,'data-preproc'];
proc.params.data_dir = [dirm,'data-preproc'];
proc.params.output_dir = [dirm,'results'];

% load subject info
% load(proc.params.qc.file);
% qc_score = cell2mat(qc_eeg(:,1));
% idx = qc_score>1; %=proc.params.qc.cutoff;
% subjects = qc_eeg(idx,2);

% get rid of this
subjects = {'2751'};

freqs = 2:0.05:45;

%% loop over subjects

for i = 1%:length(subjects)

    subject = subjects{i};
    fprintf('Analyzing subject %i of %i\n', i, length(subjects))
    
    outdir = sprintf('%s/%s/%s', preproc.params.root_dir, preproc.params.output_dir, subject);
    figdir = sprintf('%s/figures', outdir);
    
    data_file = sprintf('%s/processing_results_eeg_artfrej.mat',outdir);
    
    if ~exist(data_file, 'file')
        fprintf('No data for subject %s\n', subject);
        p = p + 1;
        continue;
    end
    
    % load data
    load(data_file);
    
    % interpolate data
    cfg             = [];
    cfg.methode     = 'linear';
    cfg.prewindow   = 0.002;
    cfg.postwindow   = 0.002; 
    data_int = ft_interpolatenan(cfg, data.eeg.ft);
    
    % find trials
    sim = get_simulation_events_eeg( outdir, preproc, proc, data );
    trials = get_trials_eeg( data, sim, params.eeg.erp );
    Ts = 1000 / data.eeg.ft.fsample;
    nchannels = size(data_int.trials{1},1);
    
    % define trials
    cfg = [];
    cfg.trl = trials.left_change.trl(:,1:3);
    data_trials_left = ft_redefinetrial(cfg,data_int);
    
    cfg = [];
    cfg.trl = trials.right_change.trl(:,1:3);
    data_trials_right = ft_redefinetrial(cfg,data_int);
    
    % run fft on full trial series
    cfg         = [];
    cfg.method  = 'pmtm'; %'pmtm';
    cfg.toi     = [0,5]; % s
    cfg.freq    = freqs;
    cfg.fcor    = false;
%     cfg.channels = 28;
    
    results_fft_left = analyze_eeg_fft(cfg, data_trials_left);
    results_fft_right = analyze_eeg_fft(cfg, data_trials_right);
    
    cfg.toi = [-5,0];
    results_fft_bl_left = analyze_eeg_fft(cfg, data_trials_left);
    results_fft_bl_right = analyze_eeg_fft(cfg, data_trials_right);
    
    % t-test against baseline    
    ntrials_left = length(results_fft_left.powspctrm);
    sp_left = sqrt((var(cat(3,results_fft_left.powspctrm{:}),[],3) + var(cat(3,results_fft_bl_left.powspctrm{:}),[],3))/2);
    tstat_left = (mean(cat(3,results_fft_left.powspctrm{:}),3) - mean(cat(3,results_fft_bl_left.powspctrm{:}),3)) ./ ...
        sp_left * sqrt(ntrials_left/2);
    
    ntrials_right = length(results_fft_right.powspctrm);
    sp_right = sqrt((var(cat(3,results_fft_right.powspctrm{:}),[],3) + var(cat(3,results_fft_bl_right.powspctrm{:}),[],3))/2);
    tstat_right = (mean(cat(3,results_fft_right.powspctrm{:}),3) - mean(cat(3,results_fft_bl_right.powspctrm{:}),3)) ./ ...
        sp_right * sqrt(ntrials_right/2);
    
    % rightlbl
    sp_right_lbl = sqrt(((ntrials_right-1)*var(cat(3,results_fft_right.powspctrm{:}),[],3) ...
        + (ntrials_left-1)*var(cat(3,results_fft_bl_left.powspctrm{:}),[],3))/ ...
        (ntrials_right + ntrials_left - 2));
    tstat_right_lbl = (mean(cat(3,results_fft_right.powspctrm{:}),3) - mean(cat(3,results_fft_bl_left.powspctrm{:}),3)) ./ ...
        (sp_right_lbl * sqrt(1/ntrials_left + 1/ntrials_right));
    
    % leftrbl
    sp_left_rbl = sqrt(((ntrials_right-1)*var(cat(3,results_fft_right.powspctrm{:}),[],3) ...
        + (ntrials_left-1)*var(cat(3,results_fft_bl_left.powspctrm{:}),[],3))/ ...
        (ntrials_right + ntrials_left - 2));
    tstat_rbl_left = (mean(cat(3,results_fft_bl_right.powspctrm{:}),3) - mean(cat(3,results_fft_left.powspctrm{:}),3)) ./ ...
        (sp_right_lbl * sqrt(1/ntrials_left + 1/ntrials_right));
    
    % plot
    figure; hold on
    plot([freqs(1),freqs(end)],[0,0],'k')
    p(1) = plot(freqs,nanmean(tstat_left,1),'b');
    p(2)= plot(freqs,nanmean(tstat_right,1),'r');
    p(3)= plot(freqs,nanmean(tstat_right_lbl,1),'r:');
    p(4)= plot(freqs,nanmean(tstat_rbl_left,1),'b:');
    xlim([freqs(1),freqs(end)]);
    ylim([-3,3]);
    ylabel('Relative power change (T-score)')
    xlabel('Frequency (Hz)')
    legend(p, {'Left-change', 'Right-change', 'Right rel. to left bl','Right baseline rel. to left'}); legend boxoff
    
    figure; hold on;
    plot(freqs,log(mean(mean(cat(3,results_fft_bl_right.powspctrm{:}),3),1)),':r');
    plot(freqs,log(mean(mean(cat(3,results_fft_right.powspctrm{:}),3),1)),'r');
    plot(freqs,log(mean(mean(cat(3,results_fft_bl_left.powspctrm{:}),3))),'b:');
    plot(freqs,log(mean(mean(cat(3,results_fft_left.powspctrm{:}),3),1)),'b');
    ylabel('Power (10log)')
    xlabel('Frequency (Hz)')
    legend({'Right baseline','Right', 'Left baseline', 'Left'}); legend boxoff
    
    
    % add events to data vars
    data_left.event = num2cell(-1*trials.left_change.trl(:,3)/data_trials_left.fsample);
    data_right.event = num2cell(-1*trials.right_change.trl(:,3)/data_trials_right.fsample);
    
    % wavelet 
    cfg         = [];
    cfg.wavelet     = 'cmor4-1';
    cfg.foi         = freqs;
    cfg.dt          = 0.010;
    cfg.events      = [1]; 
    cfg.toi        = {[-1,3]};
    
    cfg.zscore      = 0;
    cfg.FTstruct    = 1;
    cfg.saveSpectra = 0;
    cfs_left = getSpectra(cfg, data_trials_left);
    
end 
    %% 
    
    
    % define trials
    
    % run fft per events
    
    