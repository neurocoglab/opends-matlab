function [ results, summary ] = process_erp_eeg( params, data, results, summary )
%%%%%%%%%%%%%%%%%%%
% Performs timelock analyses for event-related potentials (ERPs) from EEG data, 
% given a set of trials
%
%

if ~exist(results, 'var')
   results = []; 
end

results.subject = data.subject;

% Get trial definitions
trials = get_sim_trials_eeg( params, data );
results.eeg.erp.trials = trials;

% For preprocessing
cfg = data.eeg.cfg;

% Overtake onset
cfg.trl = trials.left_change.trl;
results.eeg.erp.left_change.timelock = get_timelock( cfg );
results.eeg.erp.left_change.trl = cfg.trl;
fprintf('...done overtake onsets.\n');

% Easy
idx = find(trials.left_change.trl(:,5)==1);
cfg.trl = trials.left_change.trl(idx,1:4);
results.eeg.erp.left_change.easy.timelock = get_timelock( cfg );
results.eeg.erp.left_change.easy.trl = cfg.trl;
fprintf('...done easy overtake onsets.\n');

% Difficult
idx = find(trials.left_change.trl(:,5)>1);
cfg.trl = trials.left_change.trl(idx,1:4);
results.eeg.erp.left_change.difficult.timelock = get_timelock( cfg );
results.eeg.erp.left_change.difficult.trl = cfg.trl;
fprintf('...done difficult overtake onsets.\n');

% Positive outcome
idx = find(trials.left_change.trl(:,6)>0);
cfg.trl = trials.left_change.trl(idx,1:4);
results.eeg.erp.left_change.positive.timelock = get_timelock( cfg );
results.eeg.erp.left_change.positive.trl = cfg.trl;
fprintf('...done overtake onsets with positive outcomes.\n');

% Negative outcome
idx = find(trials.left_change.trl(:,6)<0);
cfg.trl = trials.left_change.trl(idx,1:4);
results.eeg.erp.left_change.negative.timelock = get_timelock( cfg );
results.eeg.erp.left_change.negative.trl = cfg.trl;
fprintf('...done overtake onsets with negative outcomes.\n');

% Overtake offset
cfg.trl = trials.right_change.trl;
results.eeg.erp.right_change.timelock = get_timelock( cfg );
fprintf('...done overtake offsets.\n');
results.eeg.erp.right_change.trl = cfg.trl;

% Easy
idx = find(trials.right_change.trl(:,5)==1);
cfg.trl = trials.right_change.trl(idx,1:4);
results.eeg.erp.right_change.easy.timelock = get_timelock( cfg );
results.eeg.erp.right_change.easy.trl = cfg.trl;
fprintf('...done easy overtake offsets.\n');

% Difficult
idx = find(trials.right_change.trl(:,5)>1);
cfg.trl = trials.right_change.trl(idx,1:4);
results.eeg.erp.right_change.difficult.timelock = get_timelock( cfg );
results.eeg.erp.right_change.difficult.trl = cfg.trl;
fprintf('...done difficult overtake offsets.\n');

% Positive outcome
idx = find(trials.right_change.trl(:,6)>0);
cfg.trl = trials.right_change.trl(idx,1:4);
results.eeg.erp.right_change.positive.timelock = get_timelock( cfg );
results.eeg.erp.right_change.positive.trl = cfg.trl;
fprintf('...done overtake offsets with positive outcomes.\n');

% Negative outcome
idx = find(trials.right_change.trl(:,6)<0);
cfg.trl = trials.right_change.trl(idx,1:4);
results.eeg.erp.right_change.negative.timelock = get_timelock( cfg );
results.eeg.erp.right_change.negative.trl = cfg.trl;
fprintf('...done overtake offsets with negative outcomes.\n');


% Update summary
summary = update_erp_summary( params, results, summary );


    function result = get_timelock( cfg, data )
        
        if size(cfg.trl,1) < params.eeg.erp.min_trials
            fprintf('  Not enough trials found (%d < %d)!', size(cfg.trl,1), params.eeg.erp.min_trials);
            result = [];
            return
        end
        
        % Load data into trials
        cfg.demean = 'yes';
        cfg.baselinewindow = [-.3 0];
        [~,data_erp] = evalc('ft_preprocessing(cfg);');
        
        % Select channels
%         [~,data_erp] = evalc('ft_selectdata(cfg3, data_erp);');
        if ~isempty(cfg3.badchannel)
            [~,data_erp] = evalc('ft_channelrepair(cfg3, data_erp);');
        end
        
        % Remove artifacts
        [~,data_erp] = evalc('ft_rejectartifact(cfg, data_erp);');
        
        % Interpolate over artifacts
        cfg_interp = [];
        cfg_interp.method = 'linear';
        cfg_interp.prewindow = 0.004;
        cfg_interp.postwindow = 0.004;
        cfg_interp.feedback = 'no';
        
        [~,data_erp] = evalc('ft_interpolatenan_corrected(cfg_interp, data_erp);');
        [~,data_erp] = evalc('ft_resampledata(cfg2, data_erp);');
        
        data_erp = demean_trials(data_erp, [-0.3 0.0]);

        cfg_erp = [];
        cfg_erp.trials = 'all';
        cfg_erp.vartrllength = 2;

        [~,result] = evalc('ft_timelockanalysis(cfg_erp, data_erp);');
        
%         cfg = [];
%         cfg.baseline = [-0.3 0];
%         cfg.parameter = 'trial';
% %         [~,result] = evalc('ft_timelockbaseline(cfg, data_erp);');
% %         [~,result] = evalc('ft_timelockanalysis(cfg_erp, result);');
%         [~,result] = evalc('ft_timelockbaseline(cfg, result);');
%         [~,result] = evalc('ft_timelockanalysis(cfg_erp, result);');
                
    end

    function timelock = demean_trials( timelock, baseline )
        
        tbeg = nearest(timelock.time{1}, baseline(1));
        tend = nearest(timelock.time{1}, baseline(2));
        
        N_trial = length(timelock.trial);
        N_ch = size(timelock.trial{1}, 1);
        N_t = size(timelock.trial{1}, 2);
        
        T = zeros(N_ch, N_t);
        denom = zeros(N_ch, N_t);
        
        for j = 1 : N_trial
            for k = 1 : N_ch
                timelock.trial{j}(k,:) = timelock.trial{j}(k,:) - nanmean(timelock.trial{j}(k,tbeg:tend));
                idx_keep = ~isnan(timelock.trial{j}(k,:));
                T(k,idx_keep) = T(k,idx_keep) + timelock.trial{j}(k,idx_keep);
                denom(k,idx_keep) = denom(k,idx_keep) + 1;
            end
        end
        
        T = T./denom;
        a=0;
        
    end



end