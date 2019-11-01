function [ stats ] = analyze_eeg_erp( summary, params )
%%%%%%5%
% analyze_eeg_erp - Performs statistical analyses on EEG ERP data, which
% must already have been processed with processing_eeg.m
% 
% 

stats = [];

fprintf(' Analyzing ERP results\n');

% Set up configuration
cfg = [];
cfg.method = 'montecarlo';
% cfg.statistic = 'ft_statfun_depsamplesT';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05; 
cfg.clusterstatistic = 'maxsum';
cfg.alpha = 0.025;
cfg.numrandomization = 100;
cfg.ivar = 1;
cfg.parameter = 'avg';
cfg.latency = [-0.5 1]; % 'all';
cfg.spmversion = 'spm12';
cfg_nbr = [];
% cfg_nbr.method = 'template';
% cfg_nbr.template = params.eeg.template;
cfg_nbr.method = 'distance';
cfg_nbr.layout = 'biosemi64.lay';
cfg_nbr.channel = params.eeg.channels;
[~,cfg.neighbours] = evalc('ft_prepare_neighbours(cfg_nbr);');

% Passing onset v. zero
[result, timelock_ga] = cluster_ttest_onesample(summary.erp.left_change.tlocked, cfg);
stats.left_change.timelockstats = result;
stats.left_change.timelock_avr = timelock_ga;
fprintf(' ...done passing onset v. zero.\n');

% Passing onset by difficulty
[result, timelock_ga] = cluster_ttest_twosample(summary.erp.left_change.easy.tlocked, summary.erp.left_change.difficult.tlocked, cfg);
stats.left_change.diff.levels = [{'Easy'},{'Difficult'}];
stats.left_change.diff.timelockstats = result;
stats.left_change.diff.timelock_avr = timelock_ga;
fprintf(' ...done passing onset difficulty.\n');

% Passing onset by difficulty
[result, timelock_ga] = cluster_ttest_twosample(summary.erp.left_change.positive.tlocked, summary.erp.left_change.negative.tlocked, cfg);
stats.left_change.outcome.levels = [{'Positive'},{'Negative'}];
stats.left_change.outcome.timelockstats = result;
stats.left_change.outcome.timelock_avr = timelock_ga;
fprintf(' ...done passing onset outcome.\n');

% Passing offset v. zero
[result, timelock_ga] = cluster_ttest_onesample(summary.erp.right_change.tlocked, cfg);
stats.right_change.timelockstats = result;
stats.right_change.timelock_avr = timelock_ga;
fprintf(' ...done passing offset v. zero.\n');

% Passing offset by difficulty
[result, timelock_ga] = cluster_ttest_twosample(summary.erp.right_change.easy.tlocked, summary.erp.right_change.difficult.tlocked, cfg);
stats.right_change.diff.levels = [{'Easy'},{'Difficult'}];
stats.right_change.diff.timelockstats = result;
stats.right_change.diff.timelock_avr = timelock_ga;
fprintf(' ...done passing offset difficulty.\n');

% Passing offset by difficulty
[result, timelock_ga] = cluster_ttest_twosample(summary.erp.right_change.positive.tlocked, summary.erp.right_change.negative.tlocked, cfg);
stats.right_change.outcome.levels = [{'Positive'},{'Negative'}];
stats.right_change.outcome.timelockstats = result;
stats.right_change.outcome.timelock_avr = timelock_ga;
fprintf(' ...done passing offset outcome.\n');

fprintf(' Done analyzing ERP results.\n');


    function [result, timelock_ga] = cluster_ttest_onesample( timelock, cfg )
        
        keep = true(length(timelock),1);
        for ii = 1 : length(timelock)
            keep(ii) = ~isempty(timelock{ii});
        end
        
        timelock = timelock(keep);
        N_subj = length(timelock);
        
        isubs = 1:N_subj;
        
        cfg2 = [];
        cfg2.channel = 'all';
        [~,timelock_ga] = evalc('ft_timelockgrandaverage(cfg2, timelock{:})');

        cfg.design = [isubs,isubs;ones(1,N_subj),ones(1,N_subj)*2];
%         cfg.statistic = 'ft_statfun_depsamplesT';
        cfg.uvar = 1;
        cfg.ivar = 2;
        
        timelock0 = timelock{1}; % Compare to zero
        timelock0.avg(:) = 0;
        timelock0.var(:) = 0;
        timelock0 = repmat({timelock0},N_subj,1);
        
%         result = ft_timelockstatistics(cfg, timelock{:}, timelock0{:});
        [~,result] = evalc('ft_timelockstatistics(cfg, timelock{:}, timelock0{:})');
        

    end

    function [result, timelock_ga] = cluster_ttest_twosample( timelock1, timelock2, cfg )

        keep = true(length(timelock1),1);
        for ii = 1 : length(timelock1)
            keep(ii) = ~isempty(timelock1{ii}) && ~isempty(timelock2{ii});
        end
        
        timelock1 = timelock1(keep);
        timelock2 = timelock2(keep);
        N_subj = length(timelock1);
        
        isubs = 1:N_subj;
        
        cfg2 = [];
        cfg2.channel = 'all';
        [~,timelock_ga1] = evalc('ft_timelockgrandaverage(cfg2, timelock1{:})');
        [~,timelock_ga2] = evalc('ft_timelockgrandaverage(cfg2, timelock2{:})');
        
        timelock_ga = [{timelock_ga1},{timelock_ga2}];

        cfg.design = [isubs,isubs;ones(1,N_subj),ones(1,N_subj)*2];
        cfg.uvar = 1;
        cfg.ivar = 2;

        [~,result] = evalc('ft_timelockstatistics(cfg, timelock1{:}, timelock2{:})');
%         [~,result] = evalc('ft_timelockstatistics(cfg, timelock1{:}, timelock2{:})');

%         cfg = [];
%         cfg.baseline = [-0.3 0];
%         [~,result] = evalc('ft_timelockbaseline(cfg, result)');

    end



end