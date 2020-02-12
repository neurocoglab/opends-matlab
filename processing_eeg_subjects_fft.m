
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

addsave = '_3s';

freqs_fft = 1:0.5:45;
freqs_wav = [1:10,12:2:24,28:4:44];

%% loop over subjects

[tstat_left, tstat_right, tstat_right_lbl, tstat_rbl_left,...
    tstat_pass] = deal(cell(length(subjects),1));

[fft_left_subj, fft_right_subj, fft_leftbl_subj, fft_rightbl_subj,...
    fft_pass_subj, fft_bl_subj] = deal(cell(length(subjects),1));

for i = 1%:length(subjects)
    
    subject = subjects{i};
    fprintf('Analyzing subject %i of %i\n', i, length(subjects))
    
    outdir = sprintf('%s/%s/%s', preproc.params.root_dir, preproc.params.output_dir, subject);
    figdir = sprintf('%s/figures', outdir);
    
    if exist([outdir, '\results_fft.mat'])
        load([outdir, '\results_fft.mat'])
        load([outdir, '\results_cfs.mat'])
    else
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
        
        % ----------- run fft
        cfg         = [];
        cfg.method  = 'pmtm';
        cfg.toi     = [0,3]; % s
        cfg.freq    = freqs_fft;
        cfg.fcor    = false;
        %     cfg.channels = 28;
        
        results_fft_left = analyze_eeg_fft(cfg, data_trials_left);
        results_fft_right = analyze_eeg_fft(cfg, data_trials_right);
        %     results_fft_pass = analyze_eeg_fft(cfg, data_trials_pass);
        
        cfg.toi = [-3,0];
        results_fft_bl_left = analyze_eeg_fft(cfg, data_trials_left);
        results_fft_bl_right = analyze_eeg_fft(cfg, data_trials_right);
        
        cfg.toi = [];
        results_fft_wholepass = analyze_eeg_fft(cfg, data_trials_pass);
        results_fft_wholebl = analyze_eeg_fft(cfg, data_trials_bl);
        
        
        % ------------ save data
        save([outdir, '\results_fft.mat'], 'results_fft_left', 'results_fft_right',...
            'results_fft_bl_left', 'results_fft_bl_right',...
            'results_fft_wholepass', 'results_fft_wholebl')
        
        
        % -------- compute wavelet
        cfg         = [];
        cfg.wavelet     = 'cmor4-1';
        cfg.foi         = freqs_wav;
        cfg.dt          = 0.010;
        cfg.events      = [1];
        cfg.toi         = [-1,3];
        cfg.blevents    = 1;
        
        cfg.zscore      = 0;
        cfg.FTstruct    = 1;
        cfg.saveSpectra = 0;
        cfs_left = getSpectra_baseline(cfg, data_trials_left);
        
        cfg.toi = [-4,-1];
        cfs_left_bl = getSpectra_baseline(cfg, data_trials_left);
        
        % z-score relative to baseline
        dat = abs(cat(5,cfs_left.powspctrm{:})); % channels x events x freqs x times x trials
        datbl = abs(cat(5,cfs_left_bl.powspctrm{:}));
        
        std_bl = datbl;
        std_bl(std_bl<1e-4) = NaN;
        cfsZ_left = (dat - nanmean(datbl,4)) ./ nanstd(std_bl,0,4);
        
        figure; plot(data_trials_left.trial{37}(10,:))
        figure; plot(squeeze(nanmean(dat(1,1,:,:,37),4)))
        figure; imagesc(squeeze(cfsZ_left(1,1,:,:,37)))
        
        for tr =1:64
        figure; 
        pcolor(cfs_left.time{1},freqs_wav,squeeze(nanmean(dat(1,1,:,:,37),5))); shading flat
        caxis([-6,6]); 
        end
        
        figure; 
        pcolor(cfs_left.time{1},freqs_wav,squeeze(nanmean(nanmean(dat,5),1))); 
        shading flat
        
        
        figure; 
        pcolor(cfs_left.time{1},freqs_wav,squeeze(nanmean(nanmean(cfsZ_left,5),1))); 
        shading flat
        
    end
    
    
    % ------------ store fft results
    fft_left_subj{s} = mean(cat(3,results_fft_left.powspctrm{:}),3);
    fft_right_subj{s} = mean(cat(3,results_fft_right.powspctrm{:}),3);
    
    fft_leftbl_subj{s} = mean(cat(3,results_fft_bl_left.powspctrm{:}),3);
    fft_rightbl_subj{s} = mean(cat(3,results_fft_bl_right.powspctrm{:}),3);
    
    fft_pass_subj{s} = mean(cat(3,results_fft_wholepass.powspctrm{:}),3);
    fft_bl_subj{s} = mean(cat(3,results_fft_wholebl.powspctrm{:}),3);
    
    % ------------ compute t-stats
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

%%  plot fft

figure; hold on
plot([freqsoi(1),freqsoi(end)],[0,0],'k')
p(1) = plot(freqsoi,nanmean(cat(1,tstat_left{:}),1),'b');
p(2)= plot(freqsoi,nanmean(cat(1,tstat_right{:}),1),'r');
p(3)= plot(freqsoi,nanmean(cat(1,tstat_right_lbl{:}),1),'r:');
p(4)= plot(freqsoi,nanmean(cat(1,tstat_rbl_left{:}),1),'b:');
p(5)=  plot(freqsoi,nanmean(cat(1,tstat_pass{:}),1),'g');
xlim([freqsoi(1),freqsoi(end)]);
ylim([-3,3]);
ylabel('Relative power change (T-score)')
xlabel('Frequency (Hz)')
legend(p, {'Left-change', 'Right-change', 'Right rel. to left bl',...
    'Right baseline rel. to left','Passing'}); legend boxoff

figure; hold on;
plot(freqsoi,log10(mean(cat(1,fft_bl_right),1)),':r');
plot(freqsoi,log10(mean(cat(1,fft_right),1)),'r');
plot(freqsoi,log10(mean(cat(1,fft_bl_left),1)),'b:');
plot(freqsoi,log10(mean(cat(1,fft_left),1)),'b');
plot(freqsoi,log10(mean(cat(1,fft_wholepass),1)),'g');
plot(freqsoi,log10(mean(cat(1,fft_wholebl),1)),'g:');
ylabel('Power (10log)')
xlabel('Frequency (Hz)')
legend({'Right baseline','Right', 'Left baseline', 'Left',...
    'Whole pass event','Whole baseline'}); legend boxoff

%% plot baseline & passing idx

figure; hold on
for i = 1:size(sim.epochs.idx_passing,1)
    plot(sim.epochs.idx_passing(i,1:2), [1,1],'-or')
end
for i = 1:size(sim.epochs.idx_baseline,1)
    plot(sim.epochs.idx_baseline(i,1:2), [1,1],'-ob')
end

