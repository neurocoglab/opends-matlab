function [ summary ] = analyze_hilbert_eeg( params, summary )
% ANALYZE_HILBERT_EEG Performs statistical analyses on Hilbert envelope
%   output
%   
% Note: Hilbert processing must already have been performed
%

% 0. Construct variables from summary statistics (means) across subjects
%    and grand means across channels
N_subj = length(summary.subjects);
N_channels = length(summary.channels);
N_bands = height(summary.bands);
N_cycles = max(summary.N_cycles);
stats.mean.baseline.all = zeros(N_subj,N_bands,N_channels);
stats.mean.overtake.all = zeros(N_subj,N_bands,N_channels);
stats.mean.baseline.cycles = nan(N_subj,N_bands,N_channels,N_cycles);
stats.mean.overtake.cycles = nan(N_subj,N_bands,N_channels,N_cycles);
stats.grand_mean.baseline.all = zeros(N_subj,N_bands);
stats.grand_mean.overtake.all = zeros(N_subj,N_bands);
stats.grand_mean.baseline.cycles = zeros(N_subj,N_bands,N_cycles);
stats.grand_mean.overtake.cycles = zeros(N_subj,N_bands,N_cycles);

if params.sim.epochs.difficulty.apply
    stats.mean.overtake.difficulty = zeros(N_subj,N_bands,N_channels,size(summary.epochs.stats.overtake.difficulty.mean{1},3));
    stats.grand_mean.overtake.difficulty = zeros(N_subj,N_bands,size(summary.epochs.stats.overtake.difficulty.mean{1},3));
end
if params.sim.epochs.outcomes.apply
    stats.mean.overtake.outcomes = zeros(N_subj,N_bands,N_channels,size(summary.epochs.stats.overtake.outcomes.mean{1},3));
    stats.grand_mean.overtake.outcomes = zeros(N_subj,N_bands,size(summary.epochs.stats.overtake.outcomes.mean{1},3));
end 

for i = 1 : N_subj
    stats.mean.baseline.all(i,:,:) = summary.epochs.stats.baseline.mean{i};
    stats.grand_mean.baseline.all(i,:) = mean(summary.epochs.stats.baseline.mean{i},2,'omitmissing');
    stats.mean.overtake.all(i,:,:) = summary.epochs.stats.overtake.mean{i};
    stats.grand_mean.overtake.all(i,:) = mean(summary.epochs.stats.overtake.mean{i},2,'omitmissing');
    % Cycles
    for j = 1 : summary.N_cycles(i)
        stats.mean.baseline.cycles(i,:,:,j) = summary.epochs.stats.baseline.cycles.mean{i}(:,:,j);
        stats.grand_mean.baseline.cycles(i,:,j) = mean(summary.epochs.stats.baseline.cycles.mean{i}(:,:,j),2,'omitmissing');
        stats.mean.overtake.cycles(i,:,:,j) = summary.epochs.stats.overtake.cycles.mean{i}(:,:,j);
        stats.grand_mean.overtake.cycles(i,:,j) = mean(summary.epochs.stats.overtake.cycles.mean{i}(:,:,j),2,'omitmissing');
    end

    if params.sim.epochs.difficulty.apply
        stats.mean.overtake.difficulty(i,:,:,:) = summary.epochs.stats.overtake.difficulty.mean{i};
        stats.grand_mean.overtake.difficulty(i,:,:) = mean(summary.epochs.stats.overtake.difficulty.mean{i},2,'omitmissing');
    end
    if params.sim.epochs.outcomes.apply
        stats.mean.overtake.outcomes(i,:,:,:) = summary.epochs.stats.overtake.outcomes.mean{i};
        stats.grand_mean.overtake.outcomes(i,:,:) = mean(summary.epochs.stats.overtake.outcomes.mean{i},2,'omitmissing');
    end 
end

stats.ttest = [];
stats.labels = [];

% 1. Epochs
%    Does Hilbert envelope magnitude differ between epochs? 

for bb = 1 : N_bands
    % Strategy: compute paired t-test stats for each electrode
    %           correct for FWE with FDR < 0.05
    band = summary.bands.Band{bb};
    stats.ttest.grand_mean.(band) = [];
    stats.ttest.channels.(band) = [];

    %    1.1. Passing vs. baseline

    %    1.1.1. Grand means
    Y = squeeze(stats.grand_mean.baseline.all(:,bb));
    X = squeeze(stats.grand_mean.overtake.all(:,bb));
    [~, p, ~, stat] = ttest(X,Y);
    [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p, params.eeg.hilbert.analysis.p_crit_fdr, 'dep');
    
    result = [];
    result.data = [X,Y];
    result.labels = {'Overtake','Baseline'};
    result.means = [mean(X,'omitmissing');mean(Y,'omitmissing')];
    result.pvals = p;
    result.tstats = stat.tstat;
    result.stdev = stat.sd;
    result.df = stat.df;
    result.pvals_fdr = adj_p;
    result.crit_p_fdr = crit_p;
    result.sig_fdr = h;
    result.ci_fdr = adj_ci_cvrg;
    stats.ttest.grand_mean.(band).baseline_overtakes = result;

    % By rounds/cycles (as regression)
    X = zeros(N_subj*N_cycles*2,1);
    Y = zeros(N_subj*N_cycles*2,1);
    E = zeros(N_subj*N_cycles*2,1);
    for j = 1 : N_cycles
        idx = (j-1)*N_subj*2+1:j*N_subj*2;
        X(idx) = j;
        Y(idx(1:N_subj)) = squeeze(stats.grand_mean.baseline.cycles(:,bb,j));
        E(idx(1:N_subj)) = -1;
        Y(idx(N_subj+1:end)) = squeeze(stats.grand_mean.overtake.cycles(:,bb,j));
        E(idx(N_subj+1:end)) = 1;
    end

    T = array2table([X,Y,E],"VariableNames",{'Round','MeanHilbert','Epoch'});
    idx_outlier = abs(zscore(Y)) > params.eeg.hilbert.outlier_zscore;
    lm = fitlm(T,"MeanHilbert~Round+Epoch+Round:Epoch", "Exclude", idx_outlier);
    stats.linmod.cycles.grand_mean.(band).baseline_overtakes = lm;

    %    1.1.2. Separately for channels
    Y = squeeze(stats.mean.baseline.all(:,bb,:));
    X = squeeze(stats.mean.overtake.all(:,bb,:));
    [~, p, ~, stat] = ttest(X,Y);
    [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p, params.eeg.hilbert.analysis.p_crit_fdr, 'dep');
    
    result = [];
    result.means = [mean(X,'omitmissing');mean(Y,'omitmissing')];
    result.pvals = p;
    result.tstats = stat.tstat;
    result.stdev = stat.sd;
    result.df = stat.df;
    result.pvals_fdr = adj_p;
    result.crit_p_fdr = crit_p;
    result.sig_fdr = h;
    result.ci_fdr = adj_ci_cvrg;
    stats.ttest.channels.(band).baseline_overtakes = result;
    stats.labels.baseline_overtakes = {'Overtake','Baseline'};
    
    %    1.2. High vs. low difficulty
    if params.sim.epochs.difficulty.apply
    
        stats.ttest.channels.(band).overtake = [];
        stats.ttest.channels.(band).overtake.difficulty = {};
        stats.ttest.grand_mean.(band).overtake = [];
        stats.ttest.grand_mean.(band).overtake.difficulty = {};
        if bb == 1
            stats.labels.overtake.difficulty = {};
        end
    
        % One t-test for each pair of levels
        for i = 1 : size(stats.mean.overtake.difficulty,4)-1
            X = squeeze(stats.mean.overtake.difficulty(:,bb,:,i));
            Xg = squeeze(stats.grand_mean.overtake.difficulty(:,bb,i));
            for j = i+1 : size(stats.mean.overtake.difficulty,4)
                Y = squeeze(stats.mean.overtake.difficulty(:,bb,:,j));
                Yg = squeeze(stats.grand_mean.overtake.difficulty(:,bb,j));
                [~, p, ~, stat] = ttest(Xg,Yg);
                [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p, params.eeg.hilbert.analysis.p_crit_fdr, 'dep');
                results_i = [];
                results_i.means = [mean(X,'omitmissing');mean(Y,'omitmissing')];
                results_i.pvals = p;
                results_i.tstats = stat.tstat;
                results_i.stdev = stat.sd;
                results_i.df = stat.df;
                results_i.pvals_fdr = adj_p;
                results_i.crit_p_fdr = crit_p;
                results_i.sig_fdr = h;
                results_i.ci_fdr = adj_ci_cvrg;
                stats.ttest.grand_mean.(band).overtake.difficulty = [stats.ttest.grand_mean.(band).overtake.difficulty {results_i}];

                [~, p, ~, stat] = ttest(X,Y);
                [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p,0.05,'dep');
                results_i = [];
                results_i.means = [mean(X,'omitmissing');mean(Y,'omitmissing')];
                results_i.pvals = p;
                results_i.tstats = stat.tstat;
                results_i.stdev = stat.sd;
                results_i.df = stat.df;
                results_i.pvals_fdr = adj_p;
                results_i.crit_p_fdr = crit_p;
                results_i.sig_fdr = h;
                results_i.ci_fdr = adj_ci_cvrg;
                stats.ttest.channels.(band).overtake.difficulty = [stats.ttest.channels.(band).overtake.difficulty {results_i}];
                if bb == 1
                    label = sprintf('%s_%s', params.sim.epochs.overtake.difficulty.labels{1}, ...
                                             params.sim.epochs.overtake.difficulty.labels{2});
                    stats.labels.overtake.difficulty = [stats.labels.overtake.difficulty {label}];
                end
            end
        end
    
    end
    
    %    1.3. Negative vs. positive outcome
    if params.sim.epochs.outcomes.apply
    
        stats.ttest.grand_mean.(band).overtake = [];
        stats.ttest.grand_mean.(band).overtake.outcomes = {};
        stats.ttest.channels.(band).overtake = [];
        stats.ttest.channels.(band).overtake.outcomes = {};
        if bb == 1
            stats.labels.overtake.outcomes = {};
        end
    
        % One t-test for each pair of levels
        for i = 1 : size(stats.mean.overtake.outcomes,4)-1
            X = squeeze(stats.mean.overtake.outcomes(:,bb,:,i));
            Xg = squeeze(stats.grand_mean.overtake.outcomes(:,bb,i));
            for j = i+1 : size(stats.mean.overtake.outcomes,4)
                Y = squeeze(stats.mean.overtake.outcomes(:,bb,:,j));
                Yg = squeeze(stats.grand_mean.overtake.outcomes(:,bb,j));
                [~, p, ~, stat] = ttest(Xg,Yg);
                [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p, params.eeg.hilbert.analysis.p_crit_fdr, 'dep');
                results_i = [];
                results_i.data = [Xg,Yg];
                results_i.labels = params.sim.epochs.outcomes.labels;
                results_i.means = [mean(X,'omitmissing');mean(Y,'omitmissing')];
                results_i.pvals = p;
                results_i.tstats = stat.tstat;
                results_i.stdev = stat.sd;
                results_i.df = stat.df;
                results_i.pvals_fdr = adj_p;
                results_i.crit_p_fdr = crit_p;
                results_i.sig_fdr = h;
                results_i.ci_fdr = adj_ci_cvrg;
                stats.ttest.grand_mean.(band).overtake.outcomes = [stats.ttest.grand_mean.(band).overtake.outcomes {results_i}];

                [~, p, ~, stat] = ttest(X,Y);
                [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p, params.eeg.hilbert.analysis.p_crit_fdr, 'dep');
                results_i = [];
                results_i.data = [X,Y];
                results_i.labels = params.sim.epochs.outcomes.labels;
                results_i.means = [mean(X,'omitmissing');mean(Y,'omitmissing')];
                results_i.pvals = p;
                results_i.tstats = stat.tstat;
                results_i.stdev = stat.sd;
                results_i.df = stat.df;
                results_i.pvals_fdr = adj_p;
                results_i.crit_p_fdr = crit_p;
                results_i.sig_fdr = h;
                results_i.ci_fdr = adj_ci_cvrg;
                stats.ttest.channels.(band).overtake.outcomes = [stats.ttest.channels.(band).overtake.outcomes {results_i}];
                if bb == 1
                    label = sprintf('%s_%s', params.sim.epochs.outcomes.labels{1}, ...
                                             params.sim.epochs.outcomes.labels{2});
                    stats.labels.overtake.outcomes = [stats.labels.overtake.outcomes {label}];
                end
            end
        end
    
    end 
  
end

% 2. Rounds
%    Does Hilbert envelope differ across rounds?

varnames = {'Subject','Cycle','Hilbert'};
N_cycles = params.sim.rounds.max_cycles;

stats.mean.cycles = nan(N_subj,N_bands,N_channels,N_cycles);
stats.mean.baseline.cycles = nan(N_subj,N_bands,N_channels,N_cycles);
stats.mean.overtake.cycles = nan(N_subj,N_bands,N_channels,N_cycles);
stats.grand_mean.cycles = nan(N_subj,N_bands,N_cycles);
stats.grand_mean.baseline.cycles = nan(N_subj,N_bands,N_cycles);
stats.grand_mean.overtake.cycles = nan(N_subj,N_bands,N_cycles);

% 2.1. Extract from subject data
for i = 1 : N_subj
    %    2.1.1. Overall
    idx_ci = 1:min(N_cycles, summary.N_cycles(i));
    Y = summary.epochs.stats.cycles.mean{i}(:,:,idx_ci);
    stats.mean.cycles(i,:,:,idx_ci) = Y;
    stats.grand_mean.cycles(i,:,idx_ci) = mean(Y,2,"omitmissing");
    
    %    2.1.2. Within epochs (passing/baseline)
    Y = summary.epochs.stats.baseline.cycles.mean{i}(:,:,idx_ci);
    stats.mean.baseline.cycles(i,:,:,idx_ci) = Y;
    stats.grand_mean.baseline.cycles(i,:,idx_ci) = mean(Y,2,"omitmissing");
    Y = summary.epochs.stats.overtake.cycles.mean{i}(:,:,idx_ci);
    stats.mean.overtake.cycles(i,:,:,idx_ci) = Y;
    stats.grand_mean.overtake.cycles(i,:,idx_ci) = mean(Y,2,"omitmissing");
end

% 2.2. Analyze per band - regress over rounds

for bb = 1 : N_bands
    band = summary.bands.Band{bb};
    
    stats.glm.grand_mean.cycles.(band) = [];
    stats.glm.channels.cycles.(band) = [];

    % 2.2.1. Grand means - regression over time
    % Use a mixed model to test for individual betas
    X = [];
    Y = [];
    sids = cell(N_subj,1);
    itr = 0;
    for i = 1 : N_subj
        % ii = (i-1)*N_cycles+1;
        % idx = ii:ii+N_cycles-1;
        idx = itr + (1:summary.N_cycles(i));
        sids(idx) = summary.subjects(i);
        X(idx) = 1 : summary.N_cycles(i);
        Y(idx) = squeeze(mean(summary.epochs.stats.cycles.mean{i}(bb,:,:),2,"omitnan"));
        itr = itr + summary.N_cycles(i);
    end

    % Remove outliers where |Z| > threshold
    idx = find(abs(zscore(Y)) > params.eeg.hilbert.outlier_zscore);
    stats.glm.grand_mean.cycles.(band).outliers = 0;
    if ~isempty(idx)
        X(idx) = [];
        Y(idx) = [];
        sids(idx) = [];
        stats.glm.grand_mean.cycles.(band).outliers = length(idx);
    end
    test_table = array2table([X;Y]','VariableNames',varnames(2:3));
    test_table.(varnames{1}) = sids;
    model_spec = sprintf('%s~%s+(1|%s)+(%s-1|%s)', varnames{3}, varnames{2}, varnames{1}, varnames{2}, varnames{1});
    lm = fitlme(test_table, model_spec);
    
    stats.glm.grand_mean.cycles.(band).lm = lm;
    
    M = zeros(N_cycles,1);
    S = zeros(N_cycles,1);
    for i = 1 : N_cycles
        T = test_table(test_table.(varnames{2})==i,:);
        M(i) = mean(T.(varnames{3}),"omitnan");
        S(i) = std(T.(varnames{3}),1,"omitnan");
    end

    stats.glm.grand_mean.cycles.(band).means = M;
    stats.glm.grand_mean.cycles.(band).stds = S;
    stats.glm.grand_mean.cycles.(band).all = test_table;

    % 2.2.2. Separated by channel - regression over time
    

end

% 3. Correlation with pupil?
%    Does Hilbert envelope covary with pupil diameter and 1st temporal 
%     derivative?
%    3.1. Full time series

% Strategy: Assess linear relationship across full time series, with
% subject as random term

stats.corr.eye.diam.rho = nan(N_subj,N_bands,N_channels);
stats.corr.eye.diam.pval = nan(N_subj,N_bands,N_channels);

X = [];
Y = cell(N_bands,N_channels);
S = [];

for i = 1 : N_subj
    subject = summary.subjects{i};
    outdir = sprintf( '%s/%s', params.io.output_dir, subject );
    results_file = sprintf('%s/results_hilbert_eeg.mat', outdir);
    results_eeg = load(results_file);
    t_eeg = results_eeg.results.eeg.hilbert.envelopes.time;

    results_file = sprintf('%s/results_preproc_eye.mat',outdir);
    results_eye = load(results_file);
    t_eye = results_eye.data.eye.t;

    % Synchronize overlapping time series
    % Resample t_eye to t_eeg
    t_start = max(t_eeg(1),t_eye(1));
    t_end = min(t_eeg(end),t_eye(end));
    idx_eye = t_eye >= t_start & t_eye <= t_end;
    t_eye = t_eye(idx_eye);
    idx_eeg = t_eeg >= t_eye(1) & t_eeg <= t_eye(end);
    t_eeg = t_eeg(idx_eeg);

    pd = timeseries(results_eye.data.eye.blinks.diam(idx_eye)',t_eye);
    pd = resample(pd,t_eeg);

    X = [X;pd.Data];
    S = [S;ones(pd.Length,1)*i];

    Rho = nan(N_bands,N_channels);
    Pvals = nan(N_bands,N_channels);

    for bb = 1 : N_bands
        env_i = results_eeg.results.eeg.hilbert.envelopes.zscore{bb};
        for j = 1 : N_channels
            if ~isempty(env_i{j})
                [Rho(bb,j),Pvals(bb,j)] = corr(pd.Data,env_i{j}(idx_eeg));
                Y(bb,j) = {[Y{bb,j};env_i{j}(idx_eeg)]};
            end
        end
    end

    stats.corr.eye.diam.rho(i,:,:) = Rho;
    stats.corr.eye.diam.pvals(i,:,:) = Rho;
end

% Linear models
stats.linmod.eye.diam = cell(N_bands,N_channels);

for bb = 1 : N_bands
    for j = 1 : N_channels
        % Need double conversion for fit to work
        T = table(double(X),double(Y{bb,j}),S,'VariableNames', {'PD','Hilbert', 'Subject'});
        lme = fitlme(T,'PD~Hilbert+(1|Subject)');
        stats.linmod.eye.diam(bb,j) = {lme};
    end
end


% TODO: use this?
%    3.2. Moving window (dynamic)


%    3.3. Within epochs (passing/baseline)
for i = 1 : length(subjects)

    % PD for this subject
    Xi = X(S==i);

    % Baseline/overtake epoch indices
    idx = [summary.epochs.index.baseline{i};summary.epochs.index.overtake{i}];

    Rho = nan(N_bands,N_channels,2);
    Pvals = nan(N_bands,N_channels,2);
    
    for bb = 1 : N_bands
        for j = 1 : N_channels
            % Hilbert for this band/channel
            Yibj = Y{bb,j}(S==i);
            
            for k = 1 : 2
                 [Rho(bb,j),Pvals(bb,j)] = corr(Xi(idx(:,k)),Yibj(:,k));
            end
        end
    end

end



% Return the result
summary.stats = stats;

end

