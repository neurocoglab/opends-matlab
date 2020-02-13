

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
    
    
    % -------- compute wavelet
    if exist([outdir, '\results_wav.mat'])
        load([outdir, '\results_wav.mat'])
    else
        
        % EXAMPLE: left lane changes
        data.event = data.event_left;
        results_wav_left = process_wavelet_eeg(params.wav, data);
        
        % z-score relative to baseline
        dat = abs(squeeze(results_wav_left.powspctrm(2,:,:,:,:))); % channels x events x freqs x times x trials
        datbl = abs(squeeze(results_wav_left.powspctrm(1,:,:,:,:)));
        
        std_bl = datbl;
        std_bl(std_bl<1e-4) = NaN;
        cfsZ_left = (dat - nanmean(datbl,1)) ./ nanstd(std_bl,0,1);
        
        % save wav data
        save([outdir, '\results_wav.mat'], 'results_wav_left', 'results_wav_left_bl')
        
    end
    
    % ----------- run fft
    if exist([outdir, '\results_fft.mat'])
        load([outdir, '\results_fft.mat'])
    else
        % after event (see params file for toi)
        data.event = data.event_left;
        results_fft_left = process_fft_eeg(params.fft, data);
        
        data.event = data.event_right;
        results_fft_right = process_fft_eeg(params.fft, data);
        
        % save fft data
        save([outdir, '\results_fft.mat'], 'results_fft_left', 'results_fft_right')
    end
    
    % ------------ store fft results
    fft_left_subj{s} = squeeze(nanmean(results_fft_left.powspctrm(2,:,:,:),2));
    fft_right_subj{s} = squeeze(nanmean(results_fft_right.powspctrm(2,:,:,:),2));
    
    fft_leftbl_subj{s} = squeeze(nanmean(results_fft_left.powspctrm(1,:,:,:),1));
    fft_rightbl_subj{s} = squeeze(nanmean(results_fft_right.powspctrm(1,:,:,:),1));
    
    % ------------ compute fft t-stats
    % t-test against baseline - this should be a paired t-contrast
    ntrials_left = size(results_fft_left.powspctrm,2);
    sp_left = sqrt((nanvar(results_fft_left.powspctrm(2,:,:,:),[],2) + nanvar(results_fft_left.powspctrm(1,:,:,:),[],2))/2);
    tstat_left{s} = (nanmean(results_fft_left.powspctrm(2,:,:,:),2) - nanmean(results_fft_left.powspctrm(1,:,:,:),2)) ./ ...
        sp_left * sqrt(ntrials_left/2);
    
    ntrials_right = size(results_fft_right.powspctrm,2);
    sp_right = sqrt((nanvar(results_fft_right.powspctrm(2,:,:,:),[],2) + nanvar(results_fft_right.powspctrm(1,:,:,:),[],2))/2);
    tstat_right{s} = (nanmean(results_fft_right.powspctrm(2,:,:,:),2) - nanmean(results_fft_right.powspctrm(1,:,:,:),2)) ./ ...
        sp_right * sqrt(ntrials_right/2);
    
    % t-test against other trial parts
    % rightlbl
    sp_right_lbl = sqrt(((ntrials_right-1)*nanvar(results_fft_right.powspctrm(2,:,:,:),[],2) ...
        + (ntrials_left-1)*nanvar(results_fft_left.powspctrm(1,:,:,:),[],2))/ ...
        (ntrials_right + ntrials_left - 2));
    tstat_right_lbl(s,:) = (nanmean(results_fft_right.powspctrm(2,:,:,:),2) - nanmean(results_fft_left.powspctrm(1,:,:,:),2)) ./ ...
        (sp_right_lbl * sqrt(1/ntrials_left + 1/ntrials_right));
    
    % leftrbl
    sp_left_rbl = sqrt(((ntrials_right-1)*nanvar(results_fft_left.powspctrm(2,:,:,:),[],2) ...
        + (ntrials_left-1)*nanvar(results_fft_right.powspctrm(1,:,:,:),[],2))/ ...
        (ntrials_right + ntrials_left - 2));
    tstat_rbl_left{s} = (nanmean(results_fft_left.powspctrm(2,:,:,:),2) - nanmean(results_fft_right.powspctrm(1,:,:,:),2)) ./ ...
        (sp_right_lbl * sqrt(1/ntrials_left + 1/ntrials_right));
    
end

%% save things here

save([dirm, '\results_fft_subjects.mat'], 'fft_left_subj', 'fft_right_subj',...
    'fft_leftbl_subj', 'fft_rightbl_subj', 'fft_pass_subj', 'fft_bl_subj')

save([dirm, '\results_fft_stats_subjects.mat'], 'tstat_left', 'tstat_right',...
    'tstat_right_lbl', 'tstat_rbl_left', 'tstat_pass')

% TODO save wav and wav stats here as well
