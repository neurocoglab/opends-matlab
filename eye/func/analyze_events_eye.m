function [ summary ] = analyze_events_eye( params, summary )
% Run statistical analyses on event-related eye data
%

%% Traffic decisions
if params.eye.events.traffic_decision.apply
   
    T_out = [];
    M_out = [];
    vars_m = [{'Subject'},{'ConfidenceAll'},{'ConfidenceCorrect'},{'ConfidenceWrong'}, ...
              {'Confidence1Back'},{'Confidence2Back'}, ...
              {'RtAll'},{'RtCorrect'},{'RtWrong'}, ...
              {'Rt1Back'},{'Rt2Back'}, {'PctCorrect'}, ...
              {'PctCorrect1Back'},{'PctCorrect2Back'}];
    
    % Build summary table
    for i = 1 : length(summary.subjects)
        subject = summary.subjects{i};
        results_file = sprintf('%s/%s/%s/traffic_decisions.csv', params.io.output_dir, ...
                                                                 subject, ...
                                                                 params.sim.sub_dir);
        R = readtable( results_file );
        
        % Averages
        m_conf = mean(R.Confidence);
        m_rt = mean(R.RT);
        
        % By correct/incorrect
        idx_corr = find(R.Correct==1);
        idx_ncorr = find(R.Correct==0);
        
        if isempty(idx_corr)
            n_corr = 0;
            m_conf_corr = -1;
            m_rt_corr = -1;
        else
            n_corr = length(idx_corr); 
            x = R.Confidence(idx_corr);
            m_conf_corr = mean(x(x>-1));
            x = R.RT(idx_corr);
            m_rt_corr = mean(x(x>-1));
        end
        
        if isempty(idx_ncorr)
            m_conf_ncorr = -1;
            m_rt_ncorr = -1;
        else
            x = R.Confidence(idx_ncorr);
            m_conf_ncorr = mean(x(x>-1));
            x = R.RT(idx_ncorr);
            m_rt_ncorr = mean(x(x>-1));
        end
        
        % By 1/2 back
        idx_1b = find(R.Order==2);
        idx_2b = find(R.Order==1);
        
        pct_corr_1b = -1;
        pct_corr_2b = -1;
        
        if isempty(idx_1b)
            m_conf_1b = -1;
            m_rt_1b = -1;
        else
            if n_corr>0
                pct_corr_1b = sum(R.Order==2 & R.Correct==1)/sum(R.Order==2);
            end
            
            x = R.Confidence(idx_1b);
            m_conf_1b = mean(x(x>-1));
            x = R.RT(idx_1b);
            m_rt_1b = mean(x(x>-1));
        end
        
        if isempty(idx_2b)
            m_conf_2b = -1;
            m_rt_2b = -1;
        else
            if n_corr>0
                pct_corr_2b = sum(R.Order==1 & R.Correct==1)/sum(R.Order==1);
            end
            
            x = R.Confidence(idx_2b);
            m_conf_2b = mean(x(x>-1));
            x = R.RT(idx_2b);
            m_rt_2b = mean(x(x>-1));
        end
        
        % Write result
        M = cell2table([{subject},{m_conf},{m_conf_corr},{m_conf_ncorr}, ...
                        {m_conf_1b},{m_conf_2b}, ...
                        {m_rt},{m_rt_corr},{m_rt_ncorr}, ...
                        {m_rt_1b},{m_rt_2b}, ...
                        {100*n_corr/height(R)}, ...
                        {100*pct_corr_1b},{100*pct_corr_2b}, ...
                        ],'VariableNames', vars_m);
        
        if isempty(M_out)
            M_out = M;
        else
            M_out = [M_out;M];
        end
       
        if isempty(T_out)
            T_out = R;
        else
            try
            T_out = [T_out;R];
            catch
                a=0;
            end
        end
        
    end
    
    % Write result
    result_file = sprintf('%s/traffic_decisions.csv', params.io.results_dir);
    writetable( T_out, result_file );
    
    result_file = sprintf('%s/traffic_decisions_means.csv', params.io.results_dir);
    writetable( M_out, result_file );
    
end


%% Compare left-change, overtake, and right-change to zero
%  Using cluster-based inference

if params.eye.events.left_change.apply
    [result, timelock_ga] = cluster_ttest_baseline(summary.left_change, summary.subjects, params);
    if isempty(result)
       warning('No result for Passing Onset v. Baseline!'); 
    end
    summary.stats.left_change.subjects = result.subjects;
    summary.stats.left_change.excluded_subjects = result.excluded_subjects;
    summary.stats.left_change.timelockstats = result.stats;
    summary.stats.left_change.timelock_avr = timelock_ga;
end

if params.eye.events.overtake.apply
    [result, timelock_ga] = cluster_ttest_baseline(summary.overtake, summary.subjects, params);
    if isempty(result)
       warning('No result for Overtake v. Baseline!'); 
    end
    summary.stats.overtake.subjects = result.subjects;
    summary.stats.overtake.excluded_subjects = result.excluded_subjects;
    summary.stats.overtake.timelockstats = result.stats;
    summary.stats.overtake.timelock_avr = timelock_ga;
end

if params.eye.events.right_change.apply
    [result, timelock_ga] = cluster_ttest_baseline(summary.right_change, summary.subjects, params);
    if isempty(result)
       warning('No result for Passing Offset v. Baseline!'); 
    end
    summary.stats.right_change.subjects = result.subjects;
    summary.stats.right_change.excluded_subjects = result.excluded_subjects;
    summary.stats.right_change.timelockstats = result.stats;
    summary.stats.right_change.timelock_avr = timelock_ga;
end

if params.eye.events.traffic_decision.apply
    [result, timelock_ga] = cluster_ttest_baseline(summary.traffic_decision, summary.subjects, params);
    if isempty(result)
       warning('No result for Traffic Decision v. Baseline!'); 
    end
    summary.stats.traffic_decision.subjects = result.subjects;
    summary.stats.traffic_decision.excluded_subjects = result.excluded_subjects;
    summary.stats.traffic_decision.timelockstats = result.stats;
    summary.stats.traffic_decision.timelock_avr = timelock_ga;
end


% %% Compare Easy v. Difficult
if params.eye.events.difficulty.apply
    if params.eye.events.left_change.apply
        [result, timelock_ga] = cluster_ttest_twosample(summary.left_change, summary.left_change.diffs, summary.subjects, params);
        if isempty(result)
           warning('No result for Passing Onset x Difficulty!'); 
        end
        summary.stats.left_change.diff.subjects = result.subjects;
        summary.stats.left_change.diff.excluded_subjects = result.excluded_subjects;
        summary.stats.left_change.diff.timelockstats = result.stats;
        summary.stats.left_change.diff.timelock_avr = timelock_ga;
    end

    if params.eye.events.right_change.apply
        [result, timelock_ga] = cluster_ttest_twosample(summary.right_change, summary.right_change.diffs, summary.subjects, params);
        if isempty(result)
           warning('No result for Passing Offset x Difficulty!'); 
        end
        summary.stats.right_change.diff.subjects = result.subjects;
        summary.stats.right_change.diff.excluded_subjects = result.excluded_subjects;
        summary.stats.right_change.diff.timelockstats = result.stats;
        summary.stats.right_change.diff.timelock_avr = timelock_ga;
    end

    if params.eye.events.overtake.apply
        [result, timelock_ga] = cluster_ttest_twosample(summary.overtake, summary.overtake.diffs, summary.subjects, params);
        if isempty(result)
           warning('No result for Overtake x Difficulty!'); 
        end
        summary.stats.overtake.diff.subjects = result.subjects;
        summary.stats.overtake.diff.excluded_subjects = result.excluded_subjects;
        summary.stats.overtake.diff.timelockstats = result.stats;
        summary.stats.overtake.diff.timelock_avr = timelock_ga;
    end

end

%% Compare Positive v. Negative Outcome

% Get outcome valence
if params.eye.events.outcomes.apply
    N = length(summary.left_change.outcomes);
    groups = cell(N,1);
    for i = 1 : N
        outcomesi = summary.left_change.outcomes{i};
        grpi = zeros(length(outcomesi),1);
        grpi(outcomesi<0) = 1;
        grpi(outcomesi>0) = 2;
        groups(i) = {grpi};
    end
    [result, timelock_ga] = cluster_ttest_twosample(summary.left_change, groups, summary.subjects, params);
    if isempty(result)
       warning('No result for Passing Onset x Outcome!'); 
    end
    summary.stats.left_change.outcomes.subjects = result.subjects;
    summary.stats.left_change.outcomes.excluded_subjects = result.excluded_subjects;
    summary.stats.left_change.outcomes.timelockstats = result.stats;
    summary.stats.left_change.outcomes.timelock_avr = timelock_ga;

    N = length(summary.right_change.outcomes);
    groups = cell(N,1);
    for i = 1 : N
        outcomesi = summary.right_change.outcomes{i};
        grpi = zeros(length(outcomesi),1);
        grpi(outcomesi<0) = 1;
        grpi(outcomesi>0) = 2;
        groups(i) = {grpi};
    end
    [result, timelock_ga] = cluster_ttest_twosample(summary.right_change, groups, summary.subjects, params);
    if isempty(result)
       warning('No result for Passing Offset x Outcome!'); 
    end
    summary.stats.right_change.outcomes.subjects = result.subjects;
    summary.stats.right_change.outcomes.excluded_subjects = result.excluded_subjects;
    summary.stats.right_change.outcomes.timelockstats = result.stats;
    summary.stats.right_change.outcomes.timelock_avr = timelock_ga;

    N = length(summary.overtake.outcomes);
    groups = cell(N,1);
    for i = 1 : N
        outcomesi = summary.overtake.outcomes{i};
        grpi = zeros(length(outcomesi),1);
        grpi(outcomesi<0) = 1;
        grpi(outcomesi>0) = 2;
        groups(i) = {grpi};
    end
    [result, timelock_ga] = cluster_ttest_twosample(summary.overtake, groups, summary.subjects, params);
    if isempty(result)
       warning('No result for Overtake x Outcome!'); 
    end
    summary.stats.overtake.outcomes.subjects = result.subjects;
    summary.stats.overtake.outcomes.excluded_subjects = result.excluded_subjects;
    summary.stats.overtake.outcomes.timelockstats = result.stats;
    summary.stats.overtake.outcomes.timelock_avr = timelock_ga;

end

% Compare Traffic correct vs. incorrect
if params.eye.events.traffic_decision.correct.apply
    
    N = length(summary.traffic_decision.correct);
    groups = cell(N,1);
    for i = 1 : N
        groups(i) = summary.traffic_decision.correct(i);
    end
    [result, timelock_ga] = cluster_ttest_twosample(summary.traffic_decision, groups, ...
                                                    summary.subjects, params, ...
                                                    params.eye.events.traffic_decision.min_trials);
    if isempty(result)
       warning('No result for Traffic Decision Correct v. Incorrect!'); 
    end
    summary.stats.traffic_decision.correct.subjects = result.subjects;
    summary.stats.traffic_decision.correct.excluded_subjects = result.excluded_subjects;
    summary.stats.traffic_decision.correct.timelockstats = result.stats;
    summary.stats.traffic_decision.correct.timelock_avr = timelock_ga;

end

% Compare Traffic confident vs. nonconfident
if params.eye.events.traffic_decision.confidence.apply

    N = length(summary.traffic_decision.confidence);
    groups = cell(N,1);
    for i = 1 : N
        groups(i) = summary.traffic_decision.confidence(i);
    end
    [result, timelock_ga] = cluster_ttest_twosample(summary.traffic_decision, groups, ...
                                                    summary.subjects, params, ...
                                                    params.eye.events.traffic_decision.min_trials);
    if isempty(result)
       warning('No result for Traffic Decision Confidence!'); 
    end
    summary.stats.traffic_decision.confidence.subjects = result.subjects;
    summary.stats.traffic_decision.confidence.excluded_subjects = result.excluded_subjects;
    summary.stats.traffic_decision.confidence.timelockstats = result.stats;
    summary.stats.traffic_decision.confidence.timelock_avr = timelock_ga;
    
end

% Compare Traffic 1-back vs. 2-back
if params.eye.events.traffic_decision.nback.apply

    N = length(summary.traffic_decision.nback);
    groups = cell(N,1);
    for i = 1 : N
        groups(i) = summary.traffic_decision.nback(i);
    end
    [result, timelock_ga] = cluster_ttest_twosample(summary.traffic_decision, groups, ...
                                                    summary.subjects, params, ...
                                                    params.eye.events.traffic_decision.min_trials);
    if isempty(result)
       warning('No result for Traffic Decision N-back!'); 
    end
    summary.stats.traffic_decision.nback.subjects = result.subjects;
    summary.stats.traffic_decision.nback.excluded_subjects = result.excluded_subjects;
    summary.stats.traffic_decision.nback.timelockstats = result.stats;
    summary.stats.traffic_decision.nback.timelock_avr = timelock_ga;
   
    
end




    function [result, timelock_ga] = cluster_ttest_onesample( event, subjects, params )
        
        alpha = params.eye.events.alpha/2;
        min_trials = params.eye.events.min_trials;
        
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



    function [result, timelock_ga] = cluster_ttest_baseline( event, subjects, params )
        
        alpha = params.eye.events.alpha/2;
        min_trials = params.eye.events.min_trials;
        
        % Build Fieldtrip structure
        N_subj = length(subjects);
        result.subjects = {};
        result.excluded_subjects = {};
        data = [];
        data.label = {'Event'};
        data.trial = {};
        data.time = {event.t};
        data_bl = data;
        cfg = [];
        cfg.channel = 'all';
        cfg.vartrllength = 0;
        timelock = {}; % cell(N_subj,1);
        timelock_bl = {}; % cell(N_subj,1);
        for ii = 1 : N_subj
            data.trial = {};
            data_bl.trial = {};
            trial = event.tlocked_bl{ii};
            trial_bl = event.tlocked_bl2{ii};
            keeprows = sum(isnan(trial),2) == 0;
            trial = trial(keeprows,:);
            trial_bl = trial_bl(keeprows,:);
            if size(trial,1) > min_trials
                for jj = 1 : size(trial,1)
                    data.trial(end+1) = {squeeze(trial(jj,:))};
                    data_bl.trial(end+1) = {squeeze(trial_bl(jj,:))};
                end
                data.time = repmat({event.t},1,size(trial,1));
                data_bl.time = data.time;
                [~,X] = evalc('ft_timelockanalysis(cfg,data)');
                timelock = [timelock {X}];
                [~,X] = evalc('ft_timelockanalysis(cfg,data_bl)');
                timelock_bl = [timelock_bl {X}];
                result.subjects = [result.subjects subjects(ii)];
            else
                result.excluded_subjects = [result.excluded_subjects subjects(ii)];
            end
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

        [~,result.stats] = evalc('ft_timelockstatistics(cfg, timelock{:}, timelock_bl{:})');
  
    end



    function [result, timelock_ga] = cluster_ttest_twosample( event, groups, subjects, params, min_trials )
        
        alpha = params.eye.events.alpha/2;
        if nargin < 5
            min_trials = params.eye.events.min_trials;
        end
        
        % Build Fieldtrip structure
        result.subjects = {};
        result.excluded_subjects = {};
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
            if ~isempty(trial) && sum(keeprows) >= min_trials
                for jj = 1 : size(trial,1)
                    data1.trial(end+1) = {squeeze(trial(jj,:))};
                end
                data1.time = repmat({event.t},1,size(trial,1));
                [~,X] = evalc('ft_timelockanalysis(cfg,data1)');
                idx2 = groups{ii} == 2;
                trial = event.tlocked_bl{ii}(idx2,:);
                keeprows = sum(isnan(trial),2) == 0;
                trial = trial(keeprows,:);
                if ~isempty(trial) && sum(keeprows) >= min_trials
                    timelock1(ii) = {X};
                    for jj = 1 : size(trial,1)
                        data2.trial(end+1) = {squeeze(trial(jj,:))};
                    end
                    data2.time = repmat({event.t},1,size(trial,1));
                    [~,X] = evalc('ft_timelockanalysis(cfg,data2)');
                    timelock2(ii) = {X};
                    tokeep(end+1) = ii;
                    result.subjects = [result.subjects subjects(ii)];
                else
                    result.excluded_subjects = [result.excluded_subjects subjects(ii)];
                end
            else
                result.excluded_subjects = [result.excluded_subjects subjects(ii)];
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

        [~,result.stats] = evalc('ft_timelockstatistics(cfg, timelock1{:}, timelock2{:})');
  
    end


end

