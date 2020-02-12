function [ stats ] = analyze_events( summary, params )
% Run statistical analyses on event-related data
%

%% Compare left-change, overtake, and right-change to zero
%  Using cluster-based inference
alpha = params.events.alpha/2;
min_trials = params.events.min_trials;

[result, timelock_ga] = cluster_ttest_baseline(summary.left_change, summary.subjects, alpha);
if isempty(result)
   warning('No result for Passing Onset v. Baseline!'); 
end
stats.left_change.timelockstats = result;
stats.left_change.timelock_avr = timelock_ga;
% stats.left_change.timelock_bl_avr = {timelock_bl_ga};

[result, timelock_ga] = cluster_ttest_baseline(summary.overtake, summary.subjects, alpha);
if isempty(result)
   warning('No result for Overtake v. Baseline!'); 
end
stats.overtake.timelockstats = result;
stats.overtake.timelock_avr = timelock_ga;
% stats.overtake.timelock_bl_avr = {timelock_bl_ga};

[result, timelock_ga] = cluster_ttest_baseline(summary.right_change, summary.subjects, alpha);
if isempty(result)
   warning('No result for Passing Offset v. Baseline!'); 
end
stats.right_change.timelockstats = result;
stats.right_change.timelock_avr = timelock_ga;
% stats.right_change.timelock_bl_avr = {timelock_bl_ga};


% %% Compare Easy v. Difficult

[result, timelock_ga] = cluster_ttest_twosample(summary.left_change, summary.left_change.diffs, summary.subjects, alpha, min_trials);
if isempty(result)
   warning('No result for Passing Onset x Difficulty!'); 
end
stats.left_change.diff.timelockstats = result;
stats.left_change.diff.timelock_avr = timelock_ga;

[result, timelock_ga] = cluster_ttest_twosample(summary.right_change, summary.right_change.diffs, summary.subjects, alpha, min_trials);
if isempty(result)
   warning('No result for Passing Offset x Difficulty!'); 
end
stats.right_change.diff.timelockstats = result;
stats.right_change.diff.timelock_avr = timelock_ga;

[result, timelock_ga] = cluster_ttest_twosample(summary.overtake, summary.overtake.diffs, summary.subjects, alpha, min_trials);
if isempty(result)
   warning('No result for Overtake x Difficulty!'); 
end
stats.overtake.diff.timelockstats = result;
stats.overtake.diff.timelock_avr = timelock_ga;



%% Compare Positive v. Negative Outcome

% Get outcome valence
N = length(summary.left_change.outcomes);
groups = cell(N,1);
for i = 1 : N
    outcomesi = summary.left_change.outcomes{i};
    grpi = zeros(length(outcomesi),1);
    grpi(outcomesi<0) = 1;
    grpi(outcomesi>0) = 2;
    groups(i) = {grpi};
end
[result, timelock_ga] = cluster_ttest_twosample(summary.left_change, groups, summary.subjects, alpha, min_trials);
if isempty(result)
   warning('No result for Passing Onset x Outcome!'); 
end
stats.left_change.outcomes.timelockstats = result;
stats.left_change.outcomes.timelock_avr = timelock_ga;

N = length(summary.right_change.outcomes);
groups = cell(N,1);
for i = 1 : N
    outcomesi = summary.right_change.outcomes{i};
    grpi = zeros(length(outcomesi),1);
    grpi(outcomesi<0) = 1;
    grpi(outcomesi>0) = 2;
    groups(i) = {grpi};
end
[result, timelock_ga] = cluster_ttest_twosample(summary.right_change, groups, summary.subjects, alpha, min_trials);
if isempty(result)
   warning('No result for Passing Offset x Outcome!'); 
end
stats.right_change.outcomes.timelockstats = result;
stats.right_change.outcomes.timelock_avr = timelock_ga;

N = length(summary.overtake.outcomes);
groups = cell(N,1);
for i = 1 : N
    outcomesi = summary.overtake.outcomes{i};
    grpi = zeros(length(outcomesi),1);
    grpi(outcomesi<0) = 1;
    grpi(outcomesi>0) = 2;
    groups(i) = {grpi};
end
[result, timelock_ga] = cluster_ttest_twosample(summary.overtake, groups, summary.subjects, alpha, min_trials);
if isempty(result)
   warning('No result for Overtake x Outcome!'); 
end
stats.overtake.outcomes.timelockstats = result;
stats.overtake.outcomes.timelock_avr = timelock_ga;


    function [result, timelock_ga] = cluster_ttest_onesample( event, subjects, alpha )
        
        % Build Fieldtrip structure
        N_subj = length(subjects);
        data = [];
        data.label = {'Event'};
        data.trial = {};
        data.time = {event.t};
        cfg = [];
        cfg.channel = 'all';
        cfg.vartrllength = 0;
        timelock = cell(N_subj,1);
        for ii = 1 : N_subj
            data.trial = {};
            trial = event.tlocked_bl{ii};
            keeprows = sum(isnan(trial),2) == 0;
            trial = trial(keeprows,:);
            for jj = 1 : size(trial,1)
                data.trial(end+1) = {squeeze(trial(jj,:))};
            end
            data.time = repmat({event.t},1,size(trial,1));
            [~,X] = evalc('ft_timelockanalysis(cfg,data)');
            timelock(ii) = {X};
        end

        cfg = [];
        cfg.channel = 'all';
        [~,timelock_ga] = evalc('ft_timelockgrandaverage(cfg, timelock{:})');

        timelock0 = timelock{1}; % Compare to zero
        timelock0.avg(:) = 0;
%         timelock0.var(:) = 0;
        timelock0 = repmat({timelock0},N_subj,1);
        cfg = [];
        cfg.method = 'montecarlo';
        cfg.statistic = 'ft_statfun_indepsamplesT';
        cfg.correctm = 'cluster';
        cfg.clusteralpha = 0.05; 
        cfg.clusterstatistic = 'maxsum';
        cfg.alpha = alpha;
        cfg.numrandomization = 5000;
%         cfg.design = [ones(N_subj,1);ones(N_subj,1)*2];
        subs = 1:N_subj;
        cfg.design = [subs;subs,ones(N_subj,1)';ones(N_subj,1)'*2];
        cfg.uvar=1;
        cfg.ivar = 2;
        cfg.latency = 'all';
        cfg.spmversion = 'spm12';

        [~,result] = evalc('ft_timelockstatistics(cfg, timelock{:}, timelock0{:})');
  
    end



    function [result, timelock_ga] = cluster_ttest_baseline( event, subjects, alpha )
        
        % Build Fieldtrip structure
        N_subj = length(subjects);
        data = [];
        data.label = {'Event'};
        data.trial = {};
        data.time = {event.t};
        data_bl = data;
        cfg = [];
        cfg.channel = 'all';
        cfg.vartrllength = 0;
        timelock = cell(N_subj,1);
        timelock_bl = cell(N_subj,1);
        for ii = 1 : N_subj
            data.trial = {};
            data_bl.trial = {};
            trial = event.tlocked_bl{ii};
            trial_bl = event.tlocked_bl2{ii};
            keeprows = sum(isnan(trial),2) == 0;
            trial = trial(keeprows,:);
            trial_bl = trial_bl(keeprows,:);
            for jj = 1 : size(trial,1)
                data.trial(end+1) = {squeeze(trial(jj,:))};
                data_bl.trial(end+1) = {squeeze(trial_bl(jj,:))};
            end
            data.time = repmat({event.t},1,size(trial,1));
            data_bl.time = data.time;
            [~,X] = evalc('ft_timelockanalysis(cfg,data)');
            timelock(ii) = {X};
            [~,X] = evalc('ft_timelockanalysis(cfg,data_bl)');
            timelock_bl(ii) = {X};
        end

        cfg = [];
        cfg.channel = 'all';
        [~,timelock_ga1] = evalc('ft_timelockgrandaverage(cfg, timelock{:})');
        [~,timelock_ga2] = evalc('ft_timelockgrandaverage(cfg, timelock_bl{:})');
        timelock_ga = [{timelock_ga1},{timelock_ga2}];

        cfg = [];
        cfg.method = 'montecarlo';
        cfg.statistic = 'depsamplesT';
        cfg.correctm = 'cluster';
        cfg.clusteralpha = 0.05; 
        cfg.clusterstatistic = 'maxsum';
        cfg.alpha = alpha;
        cfg.numrandomization = 1000;
        subs = 1:N_subj;
        subs = [subs,subs]';
        cfg.design = [subs,[ones(N_subj,1);ones(N_subj,1)*2]]';
        cfg.uvar = 1;
        cfg.ivar = 2;
        cfg.latency = 'all';
        cfg.spmversion = 'spm12';

        [~,result] = evalc('ft_timelockstatistics(cfg, timelock{:}, timelock_bl{:})');
  
    end



    function [result, timelock_ga] = cluster_ttest_twosample( event, groups, subjects, alpha, min_trials )
        
        if nargin < 5
           min_trials = 0; 
        end
        
        % Build Fieldtrip structure
        N_subj = length(subjects);
        data1 = [];
        data1.label = {'Event'};
        data1.trial = {};
        data1.time = {event.t};
        data2 = data1;
        cfg = [];
        cfg.channel = 'all';
        cfg.vartrllength = 0;
        timelock1 = cell(N_subj,1);
        timelock2 = cell(N_subj,1);
        
        tokeep = [];

        for ii = 1 : N_subj
            data1.trial = {};
            data2.trial = {};
            idx1 = groups{ii} == 1;
            
            trial = event.tlocked_bl{ii}(idx1,:);
            keeprows = sum(isnan(trial),2) == 0;
            trial = trial(keeprows,:);
            % Only add if both groups have at least one sample
            if ~isempty(trial) && sum(keeprows) > min_trials
                for jj = 1 : size(trial,1)
                    data1.trial(end+1) = {squeeze(trial(jj,:))};
                end
                data1.time = repmat({event.t},1,size(trial,1));
                [~,X] = evalc('ft_timelockanalysis(cfg,data1)');
                idx2 = groups{ii} == 2;
                trial = event.tlocked_bl{ii}(idx2,:);
                keeprows = sum(isnan(trial),2) == 0;
                trial = trial(keeprows,:);
                if ~isempty(trial) && sum(keeprows) > min_trials
                    timelock1(ii) = {X};
                    for jj = 1 : size(trial,1)
                        data2.trial(end+1) = {squeeze(trial(jj,:))};
                    end
                    data2.time = repmat({event.t},1,size(trial,1));
                    [~,X] = evalc('ft_timelockanalysis(cfg,data2)');
                    timelock2(ii) = {X};
                    tokeep(end+1) = ii;
                end
            end
        end

        timelock1 = timelock1(tokeep);
        timelock2 = timelock2(tokeep);
        
        if isempty(timelock1) || isempty(timelock2)
            result = [];
            timelock_ga = [];
            warning('Not enough trials to analyze events (min=%d)!', min_trials);
            return; 
        end
        
        N_kept = length(tokeep);
        timelock_ga = cell(2,1);
        
        cfg = [];
        cfg.channel = 'all';
        [~,ga] = evalc('ft_timelockgrandaverage(cfg, timelock1{:})');
        timelock_ga(1) = {ga};
        [~,ga] = evalc('ft_timelockgrandaverage(cfg, timelock2{:})');
        timelock_ga(2) = {ga};

        cfg = [];
        cfg.method = 'montecarlo';
%         cfg.statistic = 'ft_statfun_indepsamplesT';
        cfg.statistic = 'depsamplesT';
        cfg.correctm = 'cluster';
        cfg.clusteralpha = 0.05; 
        cfg.clusterstatistic = 'maxsum';
        cfg.alpha = alpha;
        cfg.numrandomization = 1000;
        subs = 1:N_kept;
        subs = [subs,subs]';
        cfg.design = [subs,[ones(N_kept,1);ones(N_kept,1)*2]]';
%         cfg.design = [subs;subs,ones(N_kept,1);ones(N_kept,1)*2];
        cfg.uvar = 1;
        cfg.ivar = 2;
        cfg.latency = 'all';
        cfg.spmversion = 'spm12';

        [~,result] = evalc('ft_timelockstatistics(cfg, timelock1{:}, timelock2{:})');
  
    end


end

