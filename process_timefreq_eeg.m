function [ results ] = process_timefreq_eeg( data, trials, params )
%%%%%%%%%%%%%%%%%%%
% Performs time/frequency analyses for event-related EEG data, 
% given a set of trials
%

% For preprocessing
cfg = data.eeg.cfg;
f = fieldnames(params.eeg.timefreq.cfg);
for i = 1 : length(f)
     cfg.(f{i}) = params.eeg.timefreq.cfg.(f{i});
end
cfg.artfctdef.eog.artifact   = data.eeg.artifacts;
cfg.artfctdef.minaccepttim   = params.eeg.artifacts.minaccepttim_timefreq;
cfg.artfctdef.reject = 'nan';
cfg.detrend = 'yes';

% Remove extreme values
Z_outlier = params.eeg.timefreq.z_max;

% For removal of bad channels
cfg3 = [];
cfg3.badchannel = data.eeg.badchannels;
cfg3.layout = 'biosemi64.lay';
cfg_nbr = [];
cfg_nbr.method = 'distance';
cfg_nbr.layout = 'biosemi64.lay';
cfg_nbr.channel = params.eeg.channels;
[~,cfg3.neighbours] = evalc('ft_prepare_neighbours(cfg_nbr);');

% Select channels
% cfg4 = [];
% cfg4.channel = data.eeg.ft.label;

cfg_bl = cfg;

% Overtake onset
cfg.trl = trials.left_change.trl;
cfg_bl.trl = trials.left_change.trl_baseline;
cfg.toi = [trials.left_change.trl(1,3); ...
    (trials.left_change.trl(1,3) + trials.left_change.trl(1,2) - trials.left_change.trl(1,1))] ...
    / data.eeg.ft.fsample;

trialdefs = {1:length(cfg.trl)};
trialdefs(end+1) = {find(trials.left_change.trl(:,5)==1)}; % Easy
trialdefs(end+1) = {find(trials.left_change.trl(:,5)>1)};  % Difficult
trialdefs(end+1) = {find(trials.left_change.trl(:,6)>0)};  % Positive
trialdefs(end+1) = {find(trials.left_change.trl(:,6)<0)};  % Negative
for r = 1 : 8
    trialdefs(end+1) = {find(trials.left_change.trl(:,7)==r)}; % Rounds
end

fprintf(' Evaluating overtake onsets...\n');
tfs = get_timefreq( cfg, cfg_bl, trialdefs );

results.eeg.timefreq.left_change.all = tfs{1};
results.eeg.timefreq.left_change.easy = tfs{2};
results.eeg.timefreq.left_change.difficult = tfs{3};
results.eeg.timefreq.left_change.positive = tfs{4};
results.eeg.timefreq.left_change.negative = tfs{5};
results.eeg.timefreq.left_change.rounds = tfs(6:end);

fprintf(' Done overtake onsets.\n');

% Overtake offset
cfg.trl = trials.right_change.trl;
cfg_bl.trl = trials.right_change.trl_baseline;
cfg.toi = [trials.right_change.trl(1,3); ...
    (trials.right_change.trl(1,3) + trials.right_change.trl(1,2) - trials.right_change.trl(1,1))] ...
    / data.eeg.ft.fsample;

trialdefs = {1:length(cfg.trl)};
trialdefs(end+1) = {find(trials.right_change.trl(:,5)==1)}; % Easy
trialdefs(end+1) = {find(trials.right_change.trl(:,5)>1)};  % Difficult
trialdefs(end+1) = {find(trials.right_change.trl(:,6)>0)};  % Positive
trialdefs(end+1) = {find(trials.right_change.trl(:,6)<0)};  % Negative
for r = 1 : 8
    trialdefs(end+1) = {find(trials.right_change.trl(:,7)==r)}; % Rounds
end

fprintf(' Evaluating overtake offsets...\n');
tfs = get_timefreq( cfg, cfg_bl, trialdefs );
results.eeg.timefreq.right_change.all = tfs{1};
results.eeg.timefreq.right_change.easy = tfs{2};
results.eeg.timefreq.right_change.difficult = tfs{3};
results.eeg.timefreq.right_change.positive = tfs{4};
results.eeg.timefreq.right_change.negative = tfs{5};
results.eeg.timefreq.right_change.rounds = tfs(6:end);
fprintf(' Done overtake offsets.\n');



    function results = get_timefreq( cfg, cfg_bl, trialdefs )
        
        % Load data into trials
        [~,mydata] = evalc('ft_preprocessing(cfg);');
        [~,bldata] = evalc('ft_preprocessing(cfg_bl);');
        
        % Interpolate over bad channels
        if ~isempty(cfg3.badchannel)
            [~,mydata] = evalc('ft_channelrepair(cfg3, mydata);');
            [~,bldata] = evalc('ft_channelrepair(cfg3, bldata);');  
        end
        
        % Remove artifacts
        [~,mydata] = evalc('ft_rejectartifact(cfg,  mydata);');
        [~,bldata] = evalc('ft_rejectartifact(cfg_bl,  bldata);');
        
        % Interpolate over artifacts
        cfg_interp = [];
        cfg_interp.method = 'linear';
        cfg_interp.prewindow = 0.004;
        cfg_interp.postwindow = 0.004;
        cfg_interp.feedback = 'no';
        [~,mydata] = evalc('ft_interpolatenan_corrected(cfg_interp, mydata);');
        [~,bldata] = evalc('ft_interpolatenan_corrected(cfg_interp, bldata);');

        mydata.event = cell(1,size(cfg.trl,1));
        bldata.event = cell(1,size(cfg.trl,1));
        for a = 1 : size(cfg.trl,1)
            mydata.event(a) = {0};
            bldata.event(a) = {0};
        end

        [~,idx_eeg] = intersect(mydata.label, params.eeg.channels); % data.eeg.eeg_channels);
        N_chan = length(idx_eeg);
        cfg_timefreq = [];
        cfg_timefreq.foi = cfg.foi;
        cfg_timefreq.wavelet = cfg.wavelet;
        cfg_timefreq.dt = cfg.dt;
        cfg_timefreq.toi = cfg.toi;
        cfg_timefreq.Zscore = cfg.Zscore;
        cfg_timefreq.heavy = cfg.heavy;
        cfg_timefreq.channels = idx_eeg;
        
        N_trial = length(mydata.trial);
        keep_trials = 1:N_trial;
        
        % Express spectra as Z-scores relative to baseline spectra
        if cfg_timefreq.heavy
            cfs = cell(N_chan,1);
            textprogressbar(-1);
            textprogressbar(sprintf('  Computing spectra for %d channels: ',length(idx_eeg)));
            for c = 1 : N_chan
                idx_c = idx_eeg(c);
                mydata_c = mydata;
                bldata_c = bldata;
                mydata_c.label = mydata.label(idx_c);
                bldata_c.label = bldata.label(idx_c);
                jidx = [];
                for jj = 1 : N_trial
                    if sum(isnan(mydata_c.trial{jj}(:))) > (0.5 * length(mydata_c.trial{jj}(:))) || ...
                       sum(isnan(bldata_c.trial{jj}(:))) > (0.5 * length(bldata_c.trial{jj}(:)))
                           
%                        fprintf(' Data for trial %d discarded.\n', jj);
                    else
                        mydata_c.trial(jj) = {mydata_c.trial{jj}(idx_c,:)};
                        bldata_c.trial(jj) = {bldata_c.trial{jj}(idx_c,:)};
                        jidx = [jidx jj];
                    end
                end
                
                keep_trials = intersect(keep_trials, jidx);
                mydata_c.time = mydata_c.time(jidx);
                mydata_c.trial = mydata_c.trial(jidx);
                mydata_c.event = mydata_c.event(jidx);
                mydata_c.trialinfo = mydata_c.trialinfo(jidx,:);
                mydata_c.sampleinfo = mydata_c.sampleinfo(jidx,:);
                
                bldata_c.time = bldata_c.time(jidx);
                bldata_c.trial = bldata_c.trial(jidx);
                bldata_c.event = bldata_c.event(jidx);
                bldata_c.trialinfo = bldata_c.trialinfo(jidx,:);
                bldata_c.sampleinfo = bldata_c.sampleinfo(jidx,:);
                
                cfg_c = cfg_timefreq;
                cfg_c.channels = 1;
                
                cfs_c = getSpectra_baseline(cfg_c, mydata_c, bldata_c);
                cfs(c) = {cfs_c};
                textprogressbar(c/N_chan*100);
            end
            textprogressbar('.');
        else
            cfs = getSpectra_baseline(cfg_timefreq, mydata, bldata);
        end
        
        results = cell(length(trialdefs),1);
        
        for tt = 1 : length(trialdefs)
            
            idx_tt = false(N_trial,1);
            idx_tt(trialdefs{tt}) = true;
            idx_tt = idx_tt(keep_trials);
            
            result = [];
            result.timelock = [];
            if sum(idx_tt) > 2
                result.timelock.avg = zeros(N_chan,length(cfs{1}.freq),length(cfs{1}.time));
                result.timelock.var = zeros(N_chan,length(cfs{1}.freq),length(cfs{1}.time));
                result.timelock.dof = zeros(N_chan,length(cfs{1}.freq),length(cfs{1}.time));
                for c = 1 : N_chan
                    % All compared to baseline
                    if cfg_timefreq.heavy
                        X_all = squeeze(cfs{c}.powspctrm(:,idx_tt,1,:,:));
                    else
                        X_all = squeeze(cfs.powspctrm(:,idx_tt,c,:,:));
                    end
                                        
                    X_all(abs(X_all) > Z_outlier) = nan;
                    
                    result.timelock.avg(c,:,:) = squeeze(nanmean(X_all,1));
                    result.timelock.var(c,:,:) = squeeze(nanvar(X_all,1));
                    result.timelock.dof(c,:,:) = squeeze(sum(~isnan(X_all),1));

                end

                result.timelock.cfg = cfg_timefreq;
                result.timelock.label = mydata.label(idx_eeg);
                result.timelock.time = cfs{1}.time;
                result.timelock.freq = cfs{1}.freq;
                result.timelock.dimord = 'chan_freq_time';
            end
            result.trl = cfg.trl(idx_tt,:);
            result.trl_baseline = cfg_bl.trl(idx_tt,:);
            results(tt) = {result};
        end
                
    end


end