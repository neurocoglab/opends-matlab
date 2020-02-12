
%% adding fieldtrip

% this should no longer be needed once event selection is done elsewhere
% TODO - change path
addpath C:\Users\marij\Documents\MATLAB\Toolboxes\fieldtrip-20180918
ft_defaults

%% params

% TODO ------------ THIS SHOULD BE MOVED ELSEWHERE
% preproc params
load processing_eeg_params.mat
preproc = load('preproc_params_hd.mat');
proc = load('processing_params.mat');

% get rid of this
dirm = 'C:\Users\marij\Documents\Projects_Driving\';
preproc.params.root_dir = [dirm,'data-preproc'];
proc.params.data_dir = [dirm,'data-preproc'];
proc.params.output_dir = [dirm,'results'];

% get rid of this
subjects = {'2751'};


%% loop over subjects

[tstat_left, tstat_right, tstat_right_lbl, tstat_rbl_left,...
    tstat_pass] = deal(cell(length(subjects),1));

[fft_left_subj, fft_right_subj, fft_leftbl_subj, fft_rightbl_subj,...
    fft_pass_subj, fft_bl_subj] = deal(cell(length(subjects),1));

for i = 1:length(subjects)
    
    subject = subjects{i};
    fprintf('Analyzing subject %i of %i\n', i, length(subjects))
    
    % TODO -------------- THIS NEEDS CHECKING
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
    
    % TODO -------------------- THIS SHOULD BE MOVED ELSEWHERE
    % Make sure there is a field .event with the time points of
    % interest in every trial. Every row should be a different event,
    % every column a trial.
    
    % interpolate data
    data_int = ft_interpolatenan(params.tinterp, data.eeg.ft);
    
    % find trials
    sim = get_simulation_events_eeg( outdir, preproc, proc, data );
    trials = get_trials_eeg( data, sim, params.eeg.erp );
    Ts = 1000 / data.eeg.ft.fsample;
    nchannels = size(data_int.trials{1},1);
    
    % define trials
    % - left
    cfg = [];
    cfg.trl = trials.left_change.trl(:,1:3);
    data_trials_left = ft_redefinetrial(cfg,data_int);
    data_trials_left.event = num2cell(zeros(length(data_trials_left.time),1));
    
    % - right
    cfg = [];
    cfg.trl = trials.right_change.trl(:,1:3);
    data_trials_right = ft_redefinetrial(cfg,data_int);
    data_trials_right.event = num2cell(zeros(length(data_trials_right.time),1));
    
    % - passing
    trl = sim.epochs.idx_passing;
    trl(:,3) = ones(size(trl,1),1);
    cfg = [];
    cfg.trl = trl;
    data_trials_pass = ft_redefinetrial(cfg,data_int);
    data_trials_pass.event = num2cell(zeros(length(data_trials_pass.time),1));
    
    % - baseline
    trl = sim.epochs.idx_baseline;
    trl(:,3) = zeros(size(trl,1),1);
    cfg = [];
    cfg.trl = trl;
    data_trials_bl = ft_redefinetrial(cfg,data_int);
    data_trials_bl.event = num2cell(zeros(length(data_trials_bl.time),1));
    
    % --------------------
    
    
    
    % -------- compute wavelet
    if exist([outdir, '\results_wav.mat'])
        load([outdir, '\results_wav.mat'])
    else
        
        % EXAMPLE: left lane changes
        results_wav_left = get_spectra_eeg(params.wav, data_trials_left);
        
        params.wav.toi = [-4,-1];
        results_wav_left_bl = get_spectra_eeg(params.wav, data_trials_left);
        
        % z-score relative to baseline
        dat = abs(cat(5,results_wav_left.powspctrm{:})); % channels x events x freqs x times x trials
        datbl = abs(cat(5,results_wav_left_bl.powspctrm{:}));
        
        std_bl = datbl;
        std_bl(std_bl<1e-4) = NaN;
        cfsZ_left = (dat - nanmean(datbl,4)) ./ nanstd(std_bl,0,4);
        
        % save wav data
        save([outdir, '\results_wav.mat'], 'results_wav_left', 'results_wav_left_bl')
        
    end
    
    % ----------- run fft
    if exist([outdir, '\results_fft.mat'])
        load([outdir, '\results_fft.mat'])
    else
        % after event (see params file for toi)
        results_fft_left = process_fft_eeg(params.fft, data_trials_left);
        results_fft_right = process_fft_eeg(params.fft, data_trials_right);
        
        % baseline
        params.fft.toi = [-3,0];
        results_fft_bl_left = process_fft_eeg(params.fft, data_trials_left);
        results_fft_bl_right = process_fft_eeg(params.fft, data_trials_right);
        
        % whole passing event
        params.fft.toi = [];
        results_fft_wholepass = process_fft_eeg(params.fft, data_trials_pass);
        results_fft_wholebl = process_fft_eeg(params.fft, data_trials_bl);
        
        % save fft data
        save([outdir, '\results_fft.mat'], 'results_fft_left', 'results_fft_right',...
            'results_fft_bl_left', 'results_fft_bl_right',...
            'results_fft_wholepass', 'results_fft_wholebl')
    end
    
    % ------------ store fft results
    fft_left_subj{s} = mean(cat(3,results_fft_left.powspctrm{:}),3);
    fft_right_subj{s} = mean(cat(3,results_fft_right.powspctrm{:}),3);
    
    fft_leftbl_subj{s} = mean(cat(3,results_fft_bl_left.powspctrm{:}),3);
    fft_rightbl_subj{s} = mean(cat(3,results_fft_bl_right.powspctrm{:}),3);
    
    fft_pass_subj{s} = mean(cat(3,results_fft_wholepass.powspctrm{:}),3);
    fft_bl_subj{s} = mean(cat(3,results_fft_wholebl.powspctrm{:}),3);
    
    % ------------ compute fft t-stats
    % t-test against baseline - this should be a paired t-contrast
    ntrials_left = length(results_fft_left.powspctrm);
    sp_left = sqrt((var(cat(3,results_fft_left.powspctrm{:}),[],3) + var(cat(3,results_fft_bl_left.powspctrm{:}),[],3))/2);
    tstat_left{s} = (mean(cat(3,results_fft_left.powspctrm{:}),3) - mean(cat(3,results_fft_bl_left.powspctrm{:}),3)) ./ ...
        sp_left * sqrt(ntrials_left/2);
    
    ntrials_right = length(results_fft_right.powspctrm);
    sp_right = sqrt((var(cat(3,results_fft_right.powspctrm{:}),[],3) + var(cat(3,results_fft_bl_right.powspctrm{:}),[],3))/2);
    tstat_right{s} = (mean(cat(3,results_fft_right.powspctrm{:}),3) - mean(cat(3,results_fft_bl_right.powspctrm{:}),3)) ./ ...
        sp_right * sqrt(ntrials_right/2);
    
    % t-test against other trial parts
    % rightlbl
    sp_right_lbl = sqrt(((ntrials_right-1)*var(cat(3,results_fft_right.powspctrm{:}),[],3) ...
        + (ntrials_left-1)*var(cat(3,results_fft_bl_left.powspctrm{:}),[],3))/ ...
        (ntrials_right + ntrials_left - 2));
    tstat_right_lbl(s,:) = (mean(cat(3,results_fft_right.powspctrm{:}),3) - mean(cat(3,results_fft_bl_left.powspctrm{:}),3)) ./ ...
        (sp_right_lbl * sqrt(1/ntrials_left + 1/ntrials_right));
    
    % leftrbl
    sp_left_rbl = sqrt(((ntrials_right-1)*var(cat(3,results_fft_right.powspctrm{:}),[],3) ...
        + (ntrials_left-1)*var(cat(3,results_fft_bl_left.powspctrm{:}),[],3))/ ...
        (ntrials_right + ntrials_left - 2));
    tstat_rbl_left{s} = (mean(cat(3,results_fft_bl_right.powspctrm{:}),3) - mean(cat(3,results_fft_left.powspctrm{:}),3)) ./ ...
        (sp_right_lbl * sqrt(1/ntrials_left + 1/ntrials_right));
    
    % passing
    ntrials_pass = length(results_fft_wholepass.powspctrm);
    ntrials_bl = length(results_fft_wholebl.powspctrm);
    sp_passbl = sqrt(((ntrials_pass-1)*var(cat(3,results_fft_wholepass.powspctrm{:}),[],3) ...
        + (ntrials_bl-1)*var(cat(3,results_fft_wholebl.powspctrm{:}),[],3))/ ...
        (ntrials_pass + ntrials_bl - 2));
    tstat_pass{s} = (mean(cat(3,results_fft_wholepass.powspctrm{:}),3) - mean(cat(3,results_fft_wholebl.powspctrm{:}),3)) ./ ...
        (sp_passbl * sqrt(1/ntrials_pass + 1/ntrials_bl));
    
end

%% save things here

save([dirm, '\results_fft_subjects.mat'], 'fft_left_subj', 'fft_right_subj',...
    'fft_leftbl_subj', 'fft_rightbl_subj', 'fft_pass_subj', 'fft_bl_subj')

save([dirm, '\results_fft_stats_subjects.mat'], 'tstat_left', 'tstat_right',...
    'tstat_right_lbl', 'tstat_rbl_left', 'tstat_pass')

% TODO save wav and wav stats here as well
