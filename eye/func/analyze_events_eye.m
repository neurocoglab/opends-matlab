function [ summary ] = analyze_events_eye( params, summary )
% Run statistical analyses on event-related eye data
%

%% Regression analyses
if ~isempty(params.eye.events.covariates.glms)

    % Read in covariates
    % subject column should be read as string
    metadata_file = sprintf('%s/%s/%s', ...
                            params.io.input_dir, ...
                            params.io.metadata_dir, ...   
                            params.general.subject_metadata_file);
    opts = detectImportOptions(metadata_file);
    opts.VariableTypes(1) = {'char'};
    metadata = readtable( metadata_file, opts );
    
    % Analyze each GLM
    for i = 1 : length(params.eye.events.covariates.glms)
        
        glm_i = params.eye.events.covariates.glms{i};

        event_names = {'left_change','right_change'};
        event_titles = {'Overtake Onset','Overtake Offset'};
        subjects = {summary.subjects,summary.subjects};
        if params.sim.events.traffic_decision.apply
            event_names(end+1) = {'traffic_decision'};
            event_titles(end+1) = {'Traffic Decision'};
            subjects = [subjects {summary.traffic_decision.subjects}];
        end

        % Run cluster analyses
        for j = 1 : length(event_names)
            
            % Get metadata + subject ids for this GLM
            has_missing = false;
            subjects_i = {};
            ivar_i = [];
            idx_i = [];
            k = 1;
            subjects_j = subjects{j};
            for s = 1 : length(subjects_j)
                subject = subjects_j{s};
                idx = find(strcmp(metadata.SubjectID, subject));
                if ~isempty(idx)
                    val = metadata.(glm_i.ivar)(idx);
                    if ~isnan(val)
                        subjects_i(k) = {subject};
                        ivar_i(k) = val;
                        idx_i(k) = s;
                        k = k + 1;
                    else
                        has_missing = true;
                    end
                end
            end

            if has_missing
                warning('GLM "%s" [%s]: covariate %s has missing values.', event_names{j}, glm_i.name, glm_i.ivar);
            end
            
            event_name = event_names{j};
            event = summary.(event_name);
            if has_missing
                event.tlocked = event.tlocked(idx_i);
                event.tlocked_bl = event.tlocked_bl(idx_i);
                event.tlocked_bl2 = event.tlocked_bl2(idx_i);
%                 event.diffs = event.diffs(idx_i);
%                 event.outcomes = event.outcomes(idx_i);
            end
            
            [result, timelock_ga] = cluster_regression( event, ivar_i, [], ...
                                                        subjects_i, params, ...
                                                        event_name );

            glm_result = {};
            glm_result.timelock_stats = result.stats;
            glm_result.timelock_ga = timelock_ga;
            glm_result.glm = glm_i;
            glm_result.event_title = event_titles{j};
            glm_result.statistic = sprintf('T statistic (%d df)', (length(subjects_i)-2));

            summary.stats.covariate_glms.(glm_i.name).(event_name) = glm_result;
            
        end
        
    end

end


%% Traffic decisions
if params.sim.events.traffic_decision.apply
   
    T_out = [];
    M_out = [];
    vars_m = [{'Subject'},{'ConfidenceAll'},{'ConfidenceCorrect'},{'ConfidenceWrong'}, ...
              {'RtAll'},{'RtCorrect'},{'RtWrong'},{'PctCorrect'}];
    
    idx_ok = true(length(summary.subjects),1);
          
    % Build summary table
    for i = 1 : length(summary.subjects)
        subject = summary.subjects{i};
        results_file = sprintf('%s/%s/%s/traffic_decisions.csv', params.io.output_dir, ...
                                                                 subject, ...
                                                                 params.sim.sub_dir);
        if ~exist( results_file, 'file' )
           warning('Subject %s has no traffic decision data', subject);
           idx_ok(i) = false;
           continue;
        end
        
        R = readtable( results_file );
        
        % Averages
        m_conf = mean(R.Confidence);
        m_rt = mean(R.RT);
        idx_corr = find(R.Correct==1);
        idx_ncorr = find(R.Correct==0);
        
        if isempty(idx_corr)
            n_corr = 0;
            m_conf_corr = -1;
            m_rt_corr = -1;
        else
            n_corr = length(idx_corr); 
            m_conf_corr = mean(R.Confidence(idx_corr));
            m_rt_corr = mean(R.RT(idx_corr));
        end
        
        if isempty(idx_ncorr)
            m_conf_ncorr = -1;
            m_rt_ncorr = -1;
        else
            m_conf_ncorr = mean(R.Confidence(idx_ncorr));
            m_rt_ncorr = mean(R.RT(idx_ncorr));
        end
        
        M = cell2table([{subject},{m_conf},{m_conf_corr},{m_conf_ncorr}, ...
                        {m_rt},{m_rt_corr},{m_rt_ncorr},{100*n_corr/height(R)}],'VariableNames', vars_m);
        
        if isempty(M_out)
            M_out = M;
        else
            M_out = [M_out;M];
        end
       
        if isempty(T_out)
            T_out = R;
        else
            T_out = [T_out;R];
        end
        
    end
    
    % Write result
    result_file = sprintf('%s/traffic_decisions.csv', params.io.results_dir);
    writetable( T_out, result_file );
    
    result_file = sprintf('%s/traffic_decisions_means.csv', params.io.results_dir);
    writetable( M_out, result_file );
    
    % Analyse versus baseline
    [result, timelock_ga] = cluster_ttest_baseline(summary.traffic_decision, ...
                                                   summary.traffic_decision.subjects, ...
                                                   params, ...
                                                   'Traffic Decision v. Baseline');
    if isempty(result)
       warning('No result for Traffic Decision v. Baseline!'); 
    else
        summary.stats.traffic_decision.subjects = result.subjects;
        summary.stats.traffic_decision.excluded_subjects = result.excluded_subjects;
        summary.stats.traffic_decision.timelock_stats = result.stats;
        summary.stats.traffic_decision.timelock_avr = timelock_ga;

    end
    
    % GLM analysis?
    
    
    % Compare correct v. incorrect
    if params.eye.events.traffic_decision.correct.apply
        params2 = params;
        params2.eye.events.min_trials=params.eye.events.traffic_decision.min_trials;
        [result, timelock_ga] = cluster_ttest_twosample(summary.traffic_decision, ...
                                                        summary.traffic_decision.correct, ...
                                                        summary.traffic_decision.subjects, ...
                                                        {'Correct','Incorrect'}, ...
                                                        params2, ...
                                                        'Traffic Correct v. Incorrect');
        if isempty(result)
           warning('No result for Traffic Decision Correct v. Incorrect!'); 
        else
            summary.stats.traffic_decision.correct.subjects = result.subjects;
            summary.stats.traffic_decision.correct.excluded_subjects = result.excluded_subjects;
            summary.stats.traffic_decision.correct.timelock_stats = result.stats;
            summary.stats.traffic_decision.correct.timelock_avr = timelock_ga;
            
            % Analyse slope & amplitude
            if params.eye.events.tlock_params.apply
                result = get_timelocked_averages(summary.traffic_decision, ...
                                             summary.traffic_decision.correct, ...
                                             summary.traffic_decision.subjects, ...
                                             {'Correct','Incorrect'});
                summary.stats.traffic_decision.correct.parameter_stats = ...
                    parameter_ttest_twosample( result, ...
                                            params.eye.events.traffic_decision, ...
                                            params.eye.events.alpha );
    
            end
        end
    end
    
    % Compare high v. low confidence
    if params.eye.events.traffic_decision.confidence.apply
        params2 = params;
        params2.eye.events.min_trials=params.eye.events.traffic_decision.min_trials;
        [result, timelock_ga] = cluster_ttest_twosample(summary.traffic_decision, ...
                                                        summary.traffic_decision.confidence, ...
                                                        summary.traffic_decision.subjects, ...
                                                        {'HighConfidence','LowConfidence'}, ...
                                                        params2, ...
                                                        'Traffic High v. Low Confidence');
        if isempty(result)
           warning('No result for Traffic Decision High v. Low Confidence!'); 
        else
            summary.stats.traffic_decision.confidence.subjects = result.subjects;
            summary.stats.traffic_decision.confidence.excluded_subjects = result.excluded_subjects;
            summary.stats.traffic_decision.confidence.timelock_stats = result.stats;
            summary.stats.traffic_decision.confidence.timelock_avr = timelock_ga;
            
            % Analyse slope & amplitude
            if params.eye.events.tlock_params.apply
                result = get_timelocked_averages(summary.traffic_decision, ...
                                                 summary.traffic_decision.confidence, ...
                                                 summary.traffic_decision.subjects, ...
                                                 {'HighConfidence','LowConfidence'});
                summary.stats.traffic_decision.confidence.parameter_stats = ...
                    parameter_ttest_twosample( result, ...
                                            params.eye.events.traffic_decision, ...
                                            params.eye.events.alpha, ...
                                            'Traffic High v. Low Confidence');
    
            end
        end

    end

    %% Compare 1- vs. 2-back
    if params.eye.events.traffic_decision.order.apply
        params2 = params;
        params2.eye.events.min_trials=params.eye.events.traffic_decision.min_trials;
        [result, timelock_ga] = cluster_ttest_twosample(summary.traffic_decision, ...
                                                        summary.traffic_decision.order, ...
                                                        summary.traffic_decision.subjects, ...
                                                        {'TwoBack','OneBack'}, ...
                                                        params2, ...
                                                        'Traffic 1 v. 2-Back');
        if isempty(result)
           warning('No result for Traffic Decision 1- vs. 2-back!'); 
        else
            summary.stats.traffic_decision.order.subjects = result.subjects;
            summary.stats.traffic_decision.order.excluded_subjects = result.excluded_subjects;
            summary.stats.traffic_decision.order.timelock_stats = result.stats;
            summary.stats.traffic_decision.order.timelock_avr = timelock_ga;
            
            % Analyse slope & amplitude
            if params.eye.events.tlock_params.apply
                result = get_timelocked_averages(summary.traffic_decision, ...
                                                 summary.traffic_decision.order, ...
                                                 summary.traffic_decision.subjects, ...
                                                 {'TwoBack','OneBack'});
                summary.stats.traffic_decision.order.parameter_stats = ...
                    parameter_ttest_twosample( result, ...
                                            params.eye.events.traffic_decision, ...
                                            params.eye.events.alpha, ...
                                            'Traffic 1 v. 2-Back');
    
            end
        end

    end
    
    
end


%% Compare left-change, overtake, and right-change to zero
%  Using cluster-based inference

[result, timelock_ga] = cluster_ttest_baseline(summary.left_change, ...
                                               summary.subjects, ...
                                               params, ...
                                               'Left Change v. Baseline');
if isempty(result)
   warning('No result for Passing Onset v. Baseline!'); 
end
summary.stats.left_change.subjects = result.subjects;
summary.stats.left_change.excluded_subjects = result.excluded_subjects;
summary.stats.left_change.timelock_stats = result.stats;
summary.stats.left_change.timelock_avr = timelock_ga;

[result, timelock_ga] = cluster_ttest_baseline(summary.overtake, ...
                                               summary.subjects, ...
                                               params, ...
                                               'Overtake v. Baseline');
if isempty(result)
   warning('No result for Overtake v. Baseline!'); 
end
summary.stats.overtake.subjects = result.subjects;
summary.stats.overtake.excluded_subjects = result.excluded_subjects;
summary.stats.overtake.timelock_stats = result.stats;
summary.stats.overtake.timelock_avr = timelock_ga;

% Analyse slope
[result, timelock_ga] = cluster_ttest_baseline(summary.right_change, ...
                                               summary.subjects, ...
                                               params, ...
                                               'Right Change v. Baseline');
if isempty(result)
   warning('No result for Passing Offset v. Baseline!'); 
end
summary.stats.right_change.subjects = result.subjects;
summary.stats.right_change.excluded_subjects = result.excluded_subjects;
summary.stats.right_change.timelock_stats = result.stats;
summary.stats.right_change.timelock_avr = timelock_ga;

% Analyse slope

% %% Compare Easy v. Difficult
if params.eye.events.difficulty.apply
    [result, timelock_ga] = cluster_ttest_twosample(summary.left_change, ...
                                                    summary.left_change.diffs, ...
                                                    summary.subjects, ...
                                                    {'Easy','Difficult'}, ...
                                                    params, ...
                                                    'Left Change Easy v. Difficult' );
    if isempty(result)
       warning('No result for Passing Onset x Difficulty!'); 
    end
    summary.stats.left_change.diff.subjects = result.subjects;
    summary.stats.left_change.diff.excluded_subjects = result.excluded_subjects;
    summary.stats.left_change.diff.timelock_stats = result.stats;
    summary.stats.left_change.diff.timelock_avr = timelock_ga;
    
    % Analyse slope & amplitude
    if params.eye.events.tlock_params.apply
        result = get_timelocked_averages(summary.left_change, ...
                                         summary.left_change.diffs, ...
                                         summary.subjects, ...
                                         {'Easy','Difficult'});
        summary.stats.left_change.diff.parameter_stats = ...
            parameter_ttest_twosample( result, ...
                                       params.eye.events.left_change, ...
                                       params.eye.events.alpha, ...
                                       'Left Change Easy v. Difficult' );

    end
    
    [result, timelock_ga] = cluster_ttest_twosample(summary.right_change, ...
                                                    summary.right_change.diffs, ...
                                                    summary.subjects, ...
                                                    {'Easy','Difficult'}, ...
                                                    params, ...
                                                    'Right Change Easy v. Difficult');
    if isempty(result)
       warning('No result for Passing Offset x Difficulty!'); 
    end
    summary.stats.right_change.diff.subjects = result.subjects;
    summary.stats.right_change.diff.excluded_subjects = result.excluded_subjects;
    summary.stats.right_change.diff.timelock_stats = result.stats;
    summary.stats.right_change.diff.timelock_avr = timelock_ga;
    
    % Analyse slope & amplitude
    if params.eye.events.tlock_params.apply
        result = get_timelocked_averages(summary.right_change, ...
                                         summary.right_change.diffs, ...
                                         summary.subjects, ...
                                         {'Easy','Difficult'});
        summary.stats.right_change.diff.parameter_stats = ...
            parameter_ttest_twosample( result, ...
                                       params.eye.events.right_change, ...
                                       params.eye.events.alpha, ...
                                       'Right Change Easy v. Difficult' );

    end
    
    % Analyse slope
    summary.stats.right_change.diff.parameter_stats = ...
        parameter_ttest_twosample( result, params.eye.events.right_change, params.eye.events.alpha );

    [result, timelock_ga] = cluster_ttest_twosample(summary.overtake, ...
                                                    summary.overtake.diffs, ...
                                                    summary.subjects, ...
                                                    {'Easy','Difficult'}, ...
                                                    params, ...
                                                    'Overtake Easy v. Difficult');
    if isempty(result)
       warning('No result for Overtake x Difficulty!'); 
    end
    summary.stats.overtake.diff.subjects = result.subjects;
    summary.stats.overtake.diff.excluded_subjects = result.excluded_subjects;
    summary.stats.overtake.diff.timelock_stats = result.stats;
    summary.stats.overtake.diff.timelock_avr = timelock_ga;
    
    % Analyse slope & amplitude
    if params.eye.events.tlock_params.apply
        result = get_timelocked_averages(summary.overtake, ...
                                         summary.overtake.diffs, ...
                                         summary.subjects, ...
                                         {'Easy','Difficult'});
        summary.stats.overtake.diff.parameter_stats = ...
            parameter_ttest_twosample( result, ...
                                       params.eye.events.overtake, ...
                                       params.eye.events.alpha, ...
                                       'Overtake Easy v. Difficult' );

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
    [result, timelock_ga] = cluster_ttest_twosample(summary.left_change, ...
                                                    groups, ...
                                                    summary.subjects, ...
                                                    {'Positive','Negative'}, ...
                                                    params, ...
                                                    'Left Change Positive v. Negative');
    if isempty(result)
       warning('No result for Passing Onset x Outcome!'); 
    end
    summary.stats.left_change.outcomes.subjects = result.subjects;
    summary.stats.left_change.outcomes.excluded_subjects = result.excluded_subjects;
    summary.stats.left_change.outcomes.timelock_stats = result.stats;
    summary.stats.left_change.outcomes.timelock_avr = timelock_ga;
    
    % Analyse slope & amplitude
    if params.eye.events.tlock_params.apply
        result = get_timelocked_averages(summary.left_change, ...
                                         groups, ...
                                         summary.subjects, ...
                                         {'Positive','Negative'});
        summary.stats.left_change.outcomes.parameter_stats = ...
            parameter_ttest_twosample( result, ...
                                       params.eye.events.left_change, ...
                                       params.eye.events.alpha, ...
                                       'Left Change Positive v. Negative' );

    end

    N = length(summary.right_change.outcomes);
    groups = cell(N,1);
    for i = 1 : N
        outcomesi = summary.right_change.outcomes{i};
        grpi = zeros(length(outcomesi),1);
        grpi(outcomesi<0) = 1;
        grpi(outcomesi>0) = 2;
        groups(i) = {grpi};
    end
    [result, timelock_ga] = cluster_ttest_twosample(summary.right_change, ...
                                                    groups, ...
                                                    summary.subjects, ...
                                                    {'Positive','Negative'}, ...
                                                    params, ...
                                                    'Right Change Positive v. Negative');
    if isempty(result)
       warning('No result for Passing Offset x Outcome!'); 
    end
    summary.stats.right_change.outcomes.subjects = result.subjects;
    summary.stats.right_change.outcomes.excluded_subjects = result.excluded_subjects;
    summary.stats.right_change.outcomes.timelock_stats = result.stats;
    summary.stats.right_change.outcomes.timelock_avr = timelock_ga;
    
    % Analyse slope & amplitude
    if params.eye.events.tlock_params.apply
        result = get_timelocked_averages(summary.right_change, ...
                                         groups, ...
                                         summary.subjects, ...
                                         {'Positive','Negative'});
        summary.stats.right_change.outcomes.parameter_stats = ...
            parameter_ttest_twosample( result, ...
                                       params.eye.events.right_change, ...
                                       params.eye.events.alpha, ...
                                       'Right Change Positive v. Negative' );

    end

    N = length(summary.overtake.outcomes);
    groups = cell(N,1);
    for i = 1 : N
        outcomesi = summary.overtake.outcomes{i};
        grpi = zeros(length(outcomesi),1);
        grpi(outcomesi<0) = 1;
        grpi(outcomesi>0) = 2;
        groups(i) = {grpi};
    end
    [result, timelock_ga] = cluster_ttest_twosample(summary.overtake, ...
                                                    groups, ...
                                                    summary.subjects, ...
                                                    {'Positive','Negative'}, ...
                                                    params, ...
                                                    'Overtake Positive v. Negative');
    if isempty(result)
       warning('No result for Overtake x Outcome!'); 
    end
    summary.stats.overtake.outcomes.subjects = result.subjects;
    summary.stats.overtake.outcomes.excluded_subjects = result.excluded_subjects;
    summary.stats.overtake.outcomes.timelock_stats = result.stats;
    summary.stats.overtake.outcomes.timelock_avr = timelock_ga;
    
    % Analyse slope & amplitude
    if params.eye.events.tlock_params.apply
        result = get_timelocked_averages(summary.overtake, ...
                                         groups, ...
                                         summary.subjects, ...
                                         {'Positive','Negative'});
        summary.stats.overtake.outcomes.parameter_stats = ...
            parameter_ttest_twosample( result, ...
                                       params.eye.events.overtake, ...
                                       params.eye.events.alpha, ...
                                       'Overtake Positive v. Negative' );

    end

end

    function [result, timelock_ga] = cluster_ttest_onesample( event, subjects, params, event_name )
        
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
            if size(trial,1) > min_trials
                for jj = 1 : size(trial,1)
                    data.trial(end+1) = {squeeze(trial(jj,:))};
                end
                data.time = repmat({event.t},1,size(trial,1));
                [~,X] = evalc('ft_timelockanalysis(cfg,data)');
                timelock(ii) = {X};
            else

            end
        end

        if N_keep < N_subj
            warning('T-test [%s]: %d subjects excluded because N_trials < min_trials (%d)', event_name, N_subj-N_keep, min_trials);
        end
        
        result.timelock_avg.trials = timelock;

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

        [~,result.stats] = evalc('ft_timelockstatistics(cfg, timelock{:}, timelock0{:})');
  
    end



    function [result, timelock_ga, timelock_subjects] = cluster_ttest_baseline( event, subjects, params, event_name )
        
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
        idx_keep = [];
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
                idx_keep = [idx_keep;ii];
            else
                result.excluded_subjects = [result.excluded_subjects subjects(ii)];
            end
        end

        N_keep = length(idx_keep);
        if N_keep < N_subj
            warning('T-test [%s]: %d subjects excluded because N_trials < min_trials (%d)', event_name, N_subj-N_keep, min_trials);
        end

        result.timelock_avg = [];
        result.timelock_avg.trials = timelock;
        result.timelock_avg.baseline = timelock_bl;
        
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
        subs = idx_keep';
        subs = [subs,subs]';
        cfg.design = [subs,[ones(N_keep,1);ones(N_keep,1)*2]]';
        cfg.uvar = 1;
        cfg.ivar = 2;
        cfg.latency = 'all';
        cfg.spmversion = 'spm12';

        [~,result.stats] = evalc('ft_timelockstatistics(cfg, timelock{:}, timelock_bl{:})');
  
    end


    % Performs cluster t-test analysis on two dependent samples derived from
    % subjects. Returns the result and a timelocked grand average
    function [result, timelock_ga] = cluster_ttest_twosample( event, groups, subjects, labels, params, event_name )
        
        alpha = params.eye.events.alpha/2;
        min_trials = params.eye.events.min_trials;
        
        % Build Fieldtrip structure
        result = [];
        result.labels = labels;
        result.subjects = {};
        result.excluded_subjects = {};
        result.groups = groups;
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

            % Only add if both groups have at least min_trials samples
            if ~isempty(trial) && sum(keeprows) >= min_trials
                for jj = 1 : size(trial,1)
                    data1.trial(end+1) = {squeeze(trial(jj,:))};
                end
                data1.time = repmat({event.t},1,size(trial,1));
                % Compute average time series for this subject
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

        N_keep = length(tokeep);
        if N_keep < N_subj
            warning('T-test2 [%s]: %d subjects excluded because N_trials < min_trials (%d)', event_name, N_subj-N_keep, min_trials);
        end

        timelock1 = timelock1(tokeep);
        timelock2 = timelock2(tokeep);
        
        if isempty(timelock1) || isempty(timelock2)
            result = [];
            timelock_ga = [];
            warning('Not enough trials to analyze events (min=%d)!', min_trials);
            return; 
        end
        
        result.timelock_avg = [];
        result.timelock_avg.(labels{1}) = timelock1;
        result.timelock_avg.(labels{2}) = timelock2;
        
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

    function result = get_timelocked_averages( event, groups, subjects, labels ) 
       
        % Build Fieldtrip structure
        result = [];
        result.labels = labels;
        result.subjects = {};
        result.excluded_subjects = {};
        result.groups = groups;
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
            % Only add if both groups have at least min_trials samples
            if ~isempty(trial) && sum(keeprows) > 0
                for jj = 1 : size(trial,1)
                    data1.trial(end+1) = {squeeze(trial(jj,:))};
                end
                data1.time = repmat({event.t},1,size(trial,1));
                % Compute average time series for this subject
                [~,X] = evalc('ft_timelockanalysis(cfg,data1)');
                idx2 = groups{ii} == 2;
                trial = event.tlocked_bl{ii}(idx2,:);
                keeprows = sum(isnan(trial),2) == 0;
                trial = trial(keeprows,:);
                if ~isempty(trial) && sum(keeprows) > 0
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

        result.timelock_avg.(labels{1}) = timelock1(tokeep);
        result.timelock_avg.(labels{2}) = timelock2(tokeep);
        
        
    end

    % Performs regression t-test analysis on time series and a covariate. Returns 
    % the result and a timelocked grand average
    function [result, timelock_ga] = cluster_regression( event, ivar, cvar, subjects, params, event_name )
        
        alpha = params.eye.events.alpha/2;
        min_trials = params.eye.events.min_trials;
        
        % Build Fieldtrip structure
        result.subjects = {};
        result.excluded_subjects = {};
        N_subj = length(subjects);
        data = [];
        data.label = {'Event'};
        data.trial = {};
        data.time = {event.t};
        cfg = [];
        cfg.channel = 'all';
        cfg.vartrllength = 0;
        timelock = cell(N_subj,1);
        
        tokeep = [];

        for ii = 1 : N_subj
            
            trial = event.tlocked_bl{ii};
            keeprows = sum(isnan(trial),2) == 0;
            trial = trial(keeprows,:);
            data.trial = {};
            
            % Only add if both groups have at least min_trials samples
            if ~isempty(trial) && sum(keeprows) >= min_trials
                for jj = 1 : size(trial,1)
                    data.trial(end+1) = {squeeze(trial(jj,:))};
                end
                data.time = repmat({event.t},1,size(trial,1));
                
                % Compute average time series for this subject
                [~,X] = evalc('ft_timelockanalysis(cfg,data)');
                timelock(ii) = {X};
                tokeep(end+1) = ii;
                    
                result.subjects = [result.subjects subjects(ii)];
            else
                result.excluded_subjects = [result.excluded_subjects subjects(ii)];
            end
        end

        timelock = timelock(tokeep);
        ivar = ivar(tokeep);
        if ~isempty(cvar)
           cvar = cvar(tokeep); 
        end
        
        if isempty(timelock)
            result = [];
            timelock_ga = [];
            warning('Regression [%s]: Not enough trials to analyze events (min=%d)!', event_name, min_trials);
            return; 
        end
        
        result.timelock_avg = timelock;
        
        cfg = [];
        cfg.channel = 'all';
        [~,timelock_ga] = evalc('ft_timelockgrandaverage(cfg, timelock{:})');

        % Perform statistical analysis
        cfg = [];
        cfg.method = 'montecarlo';
        cfg.statistic = 'ft_statfun_indepsamplesregrT';
        %cfg.statistic = 'depsamplesT';
        cfg.correctm = 'cluster';
        cfg.clusteralpha = 0.05; 
        cfg.clusterstatistic = 'maxsum';
        cfg.alpha = alpha;
        cfg.numrandomization = 1000;
       
        if ~isempty(cvar)
            cfg.design = [ivar,cvar]';
            cfg.cvar = 2;
        else
            cfg.design = ivar;
        end
        
        cfg.ivar = 1;
        cfg.latency = 'all';
        cfg.spmversion = 'spm12';

        [~,result.stats] = evalc('ft_timelockstatistics(cfg, timelock{:})');
  
    end

    

    % Performs paired t-test on slope/amplitude of averaged 
    % timelocked data versus baseline
    %
    % result - the result from cluster_ttest_twosample
    % params - parameters for the analysis
    function stats = parameter_ttest_baseline( result, params, alpha, event_name )
        
        if nargin < 3
           alpha = 0.05; 
        end
        
        timelock_1 = result.timelock_avg.trials;
        timelock_2 = result.timelock_avg.baseline;
        
        tlock_params_1 = get_tlocked_parameters( timelock_1, params );
        tlock_params_2 = get_tlocked_parameters( timelock_2, params );
        
        [stats.slope.h, stats.slope.p, stats.slope.ci, stats.slope.stats] = ...
                        ttest(tlock_params_1.slopes, tlock_params_2.slopes, 'alpha', alpha);
        
        [stats.amplitude.h, stats.amplitude.p, stats.amplitude.ci, stats.amplitude.stats] = ...
                        ttest(tlock_params_1.amplitudes, tlock_params_2.amplitudes, 'alpha', alpha);
     
    end

    % Performs paired t-test on slope/amplitude of averaged 
    % timelocked data. 
    %
    % result - the result from get_timelocked_averages
    % params - parameters for the analysis
    function stats = parameter_ttest_twosample( result, params, alpha, event_name )
        
        if nargin < 3
           alpha = 0.05; 
        end
        
        timelock_1 = result.timelock_avg.(result.labels{1});
        timelock_2 = result.timelock_avg.(result.labels{2});
        
        tlock_params_1 = get_tlocked_parameters( timelock_1, params );
        tlock_params_2 = get_tlocked_parameters( timelock_2, params );
        
        stats = [];
        stats.tlock_params = {tlock_params_1,tlock_params_2};
        
        is_ok = ~isnan(tlock_params_1.slopes);
        slopes_1 = tlock_params_1.slopes(is_ok);
        is_ok = ~isnan(tlock_params_2.slopes);
        slopes_2 = tlock_params_2.slopes(is_ok);
        
        [stats.slope.h, stats.slope.p, stats.slope.ci, stats.slope.stats] = ...
                        ttest2(slopes_1, slopes_2, 'alpha', alpha);
                    
        amplitudes_1 = tlock_params_1.amplitudes(~isnan(tlock_params_1.amplitudes));
        amplitudes_2 = tlock_params_2.amplitudes(~isnan(tlock_params_2.amplitudes));
        [stats.amplitude.h, stats.amplitude.p, stats.amplitude.ci, stats.amplitude.stats] = ...
                        ttest2(amplitudes_1, amplitudes_2, 'alpha', alpha);
     
    end

end

