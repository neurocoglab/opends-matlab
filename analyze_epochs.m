function [ stats ] = analyze_epochs( summary, subject_data )
% Performs statistical analyses on epoch data
%
%

% 1. LM for baseline/passing

PD = [];
IsBaseline = [];
Subject = {};

for i = 1 : length(summary.subjects)
    PD(end+1) = mean(summary.zscore.baseline.subjects.pupil{i});
    IsBaseline(end+1) = 1;
    Subject(end+1) = summary.subjects(i);
    PD(end+1) = mean(summary.zscore.passing.subjects.pupil{i});
    IsBaseline(end+1) = 0;
    Subject(end+1) = summary.subjects(i);
end

tbl = table(Subject', IsBaseline', double(PD)', 'VariableNames', [{'Subject'},{'IsBaseline'},{'PD'}]);
stats.passing_baseline.lm = fitlm(tbl,'PD~IsBaseline');
stats.passing_baseline.data = tbl;

% 2. LM for baseline/passing and round

Rounds = [];
PD = [];
IsBaseline = [];
Subject = {};

for i = 1 : length(summary.subjects)

    for j = 1 : length(summary.zscore.cycles.baseline.subjects.pupil{i})
        m = nanmean(summary.zscore.cycles.baseline.subjects.pupil{i}{j});
        if ~isnan(m)
            PD(end+1) = m;
            IsBaseline(end+1) = 1;
            Subject(end+1) = summary.subjects(i);
            Rounds(end+1) = j;
        end
    end
    
    for j = 1 : length(summary.zscore.cycles.passing.subjects.pupil{i})
        m = nanmean(summary.zscore.cycles.passing.subjects.pupil{i}{j});
        if ~isnan(m)
            PD(end+1) = m;
            IsBaseline(end+1) = 0;
            Subject(end+1) = summary.subjects(i);
            Rounds(end+1) = j;
        end
    end
    
end

tbl = table(Subject', IsBaseline', double(PD)', Rounds', 'VariableNames', [{'Subject'},{'IsBaseline'},{'PD'},{'Round'}]);
stats.passing_baseline_rounds.lm = fitlme(tbl,'PD~Round+IsBaseline+Round*IsBaseline+(1|Subject)+(Round-1|Subject)');
stats.passing_baseline_rounds.data = tbl;

% Look at baseline and passing separately to explain interaction
tbl_bl = tbl(tbl.IsBaseline==1,:);
tbl_pass = tbl(tbl.IsBaseline==0,:);
stats.passing_baseline_rounds.lm_round_bl = fitlme(tbl_bl,'PD~Round+(1|Subject)+(Round-1|Subject)');
stats.passing_baseline_rounds.lm_round_pass = fitlme(tbl_pass,'PD~Round+(1|Subject)+(Round-1|Subject)');

% 3. LM: Passing - Baseline PD Difference v. Performance

id_subs = subject_data.UID;
subjects = sort(intersect(id_subs, summary.subjects));
DeltaPD = [];
PD_passing = [];
PD_baseline = [];
Subject = {};
Scores = [];

X = subject_data.Game_score;

for i = 1 : length(subjects)
    subject = subjects{i};
    idx = find(strcmp(id_subs, subject));
    idx2 = find(strcmp(summary.subjects, subject));
    score = X(idx(1));
    
    mp = nanmean(summary.zscore.passing.subjects.pupil{idx2});
    mb = nanmean(summary.zscore.baseline.subjects.pupil{idx2});
    if ~isnan(mp) && ~isnan(mb)
        DeltaPD(end+1) = mp - mb;
        PD_passing(end+1) = mp;
        PD_baseline(end+1) = mb;
        Subject(end+1) = summary.subjects(i);
        Scores(end+1) = score;
    end
    
end

tbl = table(Subject', double(DeltaPD)', Scores', PD_passing', PD_baseline', 'VariableNames', [{'Subject'},{'DeltaPD'},{'Scores'},{'PD_passing'},{'PD_baseline'}]);
stats.passing_baseline_dscore.lm = fitlme(tbl,'DeltaPD~Scores+(1|Subject)+(Round-1|Subject)');
stats.passing_baseline_dscore.data = tbl;

% Score with Baseline/passing
PD = [];
IsBaseline = [];
Subject = {};
Scores = [];

for i = 1 : length(subjects)
    subject = subjects{i};
    idx = find(strcmp(id_subs, subject));
    idx2 = find(strcmp(summary.subjects, subject));
    score = X(idx(1));
    
    m = nanmean(summary.zscore.baseline.subjects.pupil{idx2});
    if ~isnan(m)
        PD(end+1) = m;
        IsBaseline(end+1) = 1;
        Subject(end+1) = summary.subjects(i);
        Scores(end+1) = score;
    end
    
    m = nanmean(summary.zscore.passing.subjects.pupil{idx2});
    if ~isnan(m)
        PD(end+1) = m;
        IsBaseline(end+1) = 0;
        Subject(end+1) = summary.subjects(i);
        Scores(end+1) = score;
    end
end

tbl = table(Subject', double(PD)', IsBaseline', Scores', 'VariableNames', [{'Subject'},{'PD'},{'IsBaseline'},{'Scores'}]);
stats.passing_baseline_score.lm = fitlm(tbl,'PD~IsBaseline+Scores+IsBaseline*Scores+(1|Subject)+(Round-1|Subject)');
stats.passing_baseline_score.data = tbl;


% 4. Difficulty and Outcome

PD = [];
Outcome = {};
Difficulty = {};
Subject = {};

diffs = [1 2];
outcomes = [-1 1];

str_diffs = [{'Easy'},{'Difficult'}];
str_outcomes = [{'Positive'},{'Negative'}];

for i = 1 : length(summary.subjects)
    pd_i = summary.zscore.passing.subjects.pupil{i};
    passing_diffs = summary.passing_diffs{i}(summary.idx_passing{i});
    passing_outcomes = summary.passing_outcomes{i}(summary.idx_passing{i});
    passing_outcomes(isnan(passing_outcomes)) = 0;
    passing_outcomes(passing_outcomes<0) = -1;
    passing_outcomes(passing_outcomes>0) = 1;
    
    % Difficulty
    for d = 1 : 2
        % Outcomes
        for o = 1 : 2
            idx = passing_diffs==diffs(d) & passing_outcomes==outcomes(o);
            PD(end+1) = mean(pd_i(idx));
            Outcome(end+1) = str_outcomes(o);
            Difficulty(end+1) = str_diffs(d);
            Subject(end+1) = summary.subjects(i);
        end
        
    end
end

tbl = table(Subject', double(PD)', Outcome', Difficulty', 'VariableNames', [{'Subject'},{'PD'},{'Outcome'},{'Difficulty'}]);
stats.passing_outcome_diff.lm = fitlm(tbl,'PD~Outcome+Difficulty+Outcome*Difficulty');
stats.passing_outcome_diff.lme = fitlme(tbl,'PD~Outcome+Difficulty+Outcome*Difficulty+(1|Subject)+(Outcome-1|Subject)');
stats.passing_outcome_diff.data = tbl;

stats.passing_outcome.mean.positive = nanmean(tbl(strcmp(tbl.Outcome,'Positive'),:).PD);
stats.passing_outcome.mean.negative = nanmean(tbl(strcmp(tbl.Outcome,'Negative'),:).PD);
stats.passing_outcome.std.positive = nanstd(tbl(strcmp(tbl.Outcome,'Positive'),:).PD);
stats.passing_outcome.std.negative = nanstd(tbl(strcmp(tbl.Outcome,'Negative'),:).PD);
stats.passing_diff.mean.easy = nanmean(tbl(strcmp(tbl.Difficulty,'Easy'),:).PD);
stats.passing_diff.mean.difficult = nanmean(tbl(strcmp(tbl.Difficulty,'Difficult'),:).PD);
stats.passing_diff.std.easy = nanstd(tbl(strcmp(tbl.Difficulty,'Easy'),:).PD);
stats.passing_diff.std.difficult = nanstd(tbl(strcmp(tbl.Difficulty,'Difficult'),:).PD);


end

