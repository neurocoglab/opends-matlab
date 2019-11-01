function [ stats ] = analyze_eeg_timefreq( summary, params, output_file )
%%%%%%5%
% analyze_eeg_timefreq - Performs statistical analyses on EEG time/freq data, which
% must already have been processed with processing_eeg.m
% 
% 

stats = [];
if nargin < 3
   output_file = []; 
end

fprintf(' Analyzing time/frequency results\n');

% Set up configuration
cfg = [];
cfg.method = 'montecarlo';
% cfg.statistic = 'ft_statfun_depsamplesT';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05; 
cfg.clusterstatistic = 'maxsum';
cfg.alpha = 0.025;
cfg.numrandomization = 150;
cfg.ivar = 1;
cfg.parameter = 'avg';
cfg.latency = 'all';
cfg.spmversion = 'spm12';
cfg_nbr = [];
% cfg_nbr.method = 'template';
% cfg_nbr.template = params.eeg.template;
cfg_nbr.method = 'distance';
cfg_nbr.layout = 'biosemi64.lay';
cfg_nbr.channel = params.eeg.channels;
[~,cfg.neighbours] = evalc('ft_prepare_neighbours(cfg_nbr);');

% Passing onset v. zero

if ~isempty(output_file) && exist(sprintf('%s_left_change_all.mat', output_file), 'file')
    load(sprintf('%s_left_change_all.mat', output_file));
    stats.left_change = T;
    clear T;
    fprintf(' ...skipping passing onset v. zero.\n');
else
    [result, timelock_ga] = cluster_ttest_onesample(summary.timefreq.left_change.all.tlocked, cfg);
    stats.left_change.timelockstats = result;
    stats.left_change.timelock_avr = timelock_ga;
    if ~isempty(output_file)
        T = stats.left_change;
        save(sprintf('%s_left_change_all.mat', output_file), 'T');
        clear T;
    end
    fprintf(' ...done passing onset v. zero.\n');
end

% Passing onset by difficulty
if ~isempty(output_file) && exist(sprintf('%s_left_change_diff.mat', output_file), 'file')
    load(sprintf('%s_left_change_diff.mat', output_file));
    stats.left_change.diff = T;
    clear T;
    fprintf(' ...skipping passing onset difficulty.\n');
else
    [result, timelock_ga] = cluster_ttest_twosample(summary.timefreq.left_change.easy.tlocked, summary.timefreq.left_change.difficult.tlocked, cfg);
    stats.left_change.diff.levels = [{'Easy'},{'Difficult'}];
    stats.left_change.diff.timelockstats = result;
    stats.left_change.diff.timelock_avr = timelock_ga;
    if ~isempty(output_file)
        T = stats.left_change.diff;
        save(sprintf('%s_left_change_diff.mat', output_file), 'T');
        clear T;
    end
    fprintf(' ...done passing onset difficulty.\n');
end

% Passing onset by difficulty
if ~isempty(output_file) && exist(sprintf('%s_left_change_outcome.mat', output_file), 'file')
    load(sprintf('%s_left_change_outcome.mat', output_file));
    stats.left_change.outcome = T;
    clear T;
    fprintf(' ...skipping passing onset outcome.\n');
    
else
    [result, timelock_ga] = cluster_ttest_twosample(summary.timefreq.left_change.positive.tlocked, summary.timefreq.left_change.negative.tlocked, cfg);
    stats.left_change.outcome.levels = [{'Positive'},{'Negative'}];
    stats.left_change.outcome.timelockstats = result;
    stats.left_change.outcome.timelock_avr = timelock_ga;
    if ~isempty(output_file)
        T = stats.left_change.outcome;
        save(sprintf('%s_left_change_outcome.mat', output_file), 'T');
        clear T;
    end
    fprintf(' ...done passing onset outcome.\n');
end

% Passing offset v. zero
if ~isempty(output_file) && exist(sprintf('%s_right_change_all.mat', output_file), 'file')
    load(sprintf('%s_right_change_all.mat', output_file));
    stats.right_change = T;
    clear T;
    fprintf(' ...skipping passing offset v. zero.\n');
else
    [result, timelock_ga] = cluster_ttest_onesample(summary.timefreq.right_change.all.tlocked, cfg);
    stats.right_change.timelockstats = result;
    stats.right_change.timelock_avr = timelock_ga;
    fprintf(' ...done passing offset v. zero.\n');
    if ~isempty(output_file)
        T = stats.right_change;
        save(sprintf('%s_right_change_all.mat', output_file), 'T');
        clear T;
    end
end

% Passing offset by difficulty
if ~isempty(output_file) && exist(sprintf('%s_right_change_diff.mat', output_file), 'file')
    load(sprintf('%s_right_change_diff.mat', output_file));
    stats.right_change.diff = T;
    clear T;
    fprintf(' ...skipping passing offset difficulty.\n');
else
    [result, timelock_ga] = cluster_ttest_twosample(summary.timefreq.right_change.easy.tlocked, summary.timefreq.right_change.difficult.tlocked, cfg);
    stats.right_change.diff.levels = [{'Easy'},{'Difficult'}];
    stats.right_change.diff.timelockstats = result;
    stats.right_change.diff.timelock_avr = timelock_ga;
    if ~isempty(output_file)
        T = stats.right_change.diff;
        save(sprintf('%s_right_change_diff.mat', output_file), 'T');
        clear T;
    end
    fprintf(' ...done passing offset difficulty.\n');
end

% Passing offset by difficulty
if ~isempty(output_file) && exist(sprintf('%s_right_change_outcome.mat', output_file), 'file')
    load(sprintf('%s_right_change_outcome.mat', output_file));
    stats.right_change.outcome = T;
    clear T;
    fprintf(' ...skipping passing offset outcome.\n');
else
    [result, timelock_ga] = cluster_ttest_twosample(summary.timefreq.right_change.positive.tlocked, summary.timefreq.right_change.negative.tlocked, cfg);
    stats.right_change.outcome.levels = [{'Positive'},{'Negative'}];
    stats.right_change.outcome.timelockstats = result;
    stats.right_change.outcome.timelock_avr = timelock_ga;
    if ~isempty(output_file)
        T = stats.right_change.outcome;
        save(sprintf('%s_right_change_outcome.mat', output_file), 'T');
        clear T;
    end
    fprintf(' ...done passing offset outcome.\n');
end

fprintf(' Done analyzing time/frequency results.\n');


    function [result, timelock_ga] = cluster_ttest_onesample( timelock, cfg )
        
        keep = true(length(timelock),1);
        for ii = 1 : length(timelock)
            keep(ii) = ~isempty(timelock{ii});
            keep(ii) = keep(ii) && ~any(isinf(timelock{ii}.avg(:)));
        end
        
        timelock = timelock(keep);
        N_subj = length(timelock);
        
        isubs = 1:N_subj;
        
        cfg2 = [];
        cfg2.channel = 'all';
        cfg2.toilim = 'all';
        cfg2.foilim = 'all';
        cfg2.parameter = 'avg';
        timelock_ga = ft_freqgrandaverage(cfg2, timelock{:});
%         [~,timelock_ga] = evalc('ft_timelockgrandaverage(cfg2, timelock{:})');

        cfg.design = [isubs,isubs;ones(1,N_subj),ones(1,N_subj)*2];
%         cfg.statistic = 'ft_statfun_depsamplesT';
        cfg.uvar = 1;
        cfg.ivar = 2;
        
        timelock0 = timelock{1}; % Compare to zero
        timelock0.avg(:) = 0;
        timelock0.var(:) = 0;
        timelock0 = repmat({timelock0},N_subj,1);
        
        result = ft_freqstatistics(cfg, timelock{:}, timelock0{:});
%         [~,result] = evalc('ft_timelockstatistics(cfg, timelock{:}, timelock0{:})');
        

    end

    function [result, timelock_ga] = cluster_ttest_twosample( timelock1, timelock2, cfg )

        keep = true(length(timelock1),1);
        for ii = 1 : length(timelock1)
            keep(ii) = ~isempty(timelock1{ii}) && ~isempty(timelock2{ii});
            keep(ii) = keep(ii) && ~any(isinf(timelock1{ii}.avg(:))) && ~any(isinf(timelock2{ii}.avg(:)));
        end
        
        timelock1 = timelock1(keep);
        timelock2 = timelock2(keep);
        N_subj = length(timelock1);
        
        isubs = 1:N_subj;
        
        cfg2 = [];
        cfg2.channel = 'all';
        cfg2.toilim = 'all';
        cfg2.foilim = 'all';
        cfg2.parameter = 'avg';
        timelock_ga1 = ft_freqgrandaverage(cfg2, timelock1{:});
        timelock_ga2 = ft_freqgrandaverage(cfg2, timelock2{:});
        
%         cfg2 = [];
%         cfg2.channel = 'all';
%         [~,timelock_ga1] = evalc('ft_timelockgrandaverage(cfg2, timelock1{:})');
%         [~,timelock_ga2] = evalc('ft_timelockgrandaverage(cfg2, timelock2{:})');
        
        timelock_ga = [{timelock_ga1},{timelock_ga2}];

        cfg.design = [isubs,isubs;ones(1,N_subj),ones(1,N_subj)*2];
        cfg.uvar = 1;
        cfg.ivar = 2;

        result = ft_freqstatistics(cfg, timelock1{:}, timelock2{:});
%         [~,result] = evalc('ft_timelockstatistics(cfg, timelock1{:}, timelock2{:})');
%         [~,result] = evalc('ft_timelockstatistics(cfg, timelock1{:}, timelock2{:})');
        

    end



end