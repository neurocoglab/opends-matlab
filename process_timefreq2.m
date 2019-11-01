%%
% Get params
load processing_eeg_params2;
preproc = load('preproc_params_hd.mat');
proc = load('processing_params.mat');

addpath('../lib_areid');

load(proc.params.qc.file);
qc_score = cell2mat(qc_eeg(:,1));
idx = qc_score>1; %=proc.params.qc.cutoff;
subjects = qc_eeg(idx,2);

N_subj = length(subjects);
figure_dir = ['/Users/lpzatr/OneDrive/OneDrive UoN/OneDrive - The University' ...
             ' of Nottingham/synched/projects/driving/figures/timefreq'];

if ~exist(figure_dir, 'dir')
   mkdir(figure_dir); 
end
         
channels = [{'Fz'},{'F2'},{'AF7'},{'AF8'},{'Cz'},{'Pz'},{'FC1'},{'C3'},{'C4'}];

show_figures = 'off';

%% Read in data
results = {};

for s = 1 : N_subj
    subject = subjects{s};
    subj_dir = sprintf('%s/%s/%s', preproc.params.root_dir, preproc.params.output_dir, subject);
    
    fprintf('Reading time/freq results for %s\n', subject);
    
    results(s) = {load(sprintf('%s/processing_results_eeg_timefreq.mat',subj_dir), 'data_timefreq')};

end

fprintf('\nDone loading time/freq results.\n');


%% Extract effects matrix

channels_all = {};

timefreq = results{1}.data_timefreq.left_change{1};
foi = timefreq.cfg.foi;
toi = timefreq.cfg.toi(1) : timefreq.cfg.dt: timefreq.cfg.toi(2);

C_all = nan(length(channels), N_subj, length(foi), length(toi));
C_diff = nan(length(channels), N_subj, length(foi), length(toi));

for s = 1 : N_subj
    channels_s = results{s}.data_timefreq.left_change{1}.channels;
    channels_all = union(channels_all, channels_s);
    for c = 1 : length(channels)
        idx = find(strcmp(channels_s, channels{c}));
        if isempty(idx)
            warning('Subject %s has no channel %s', subjects{s}, channels{c});
        else
            timefreq = results{s}.data_timefreq.left_change{idx};
            C_all(c,s,:,:) = timefreq.stats.all_baseline.c;
            C_diff(c,s,:,:) = timefreq.stats.easy_diff.c;
%             
%             T_bl(c,s,:,:) = timefreq.stats.all_baseline.t;
%             T_diff(c,s,:,:) = timefreq.stats.easy_diff.t;
        end
    end
    
end

    %% Perform 2nd-level analysis
    
for c = 1 : length(channels)
        
    Cc_all = squeeze(C_all(c,:,:,:));
    Cc_diff = squeeze(C_diff(c,:,:,:));
    
    M_all = squeeze(nanmean(Cc_all,1));
    M_diff = squeeze(nanmean(Cc_diff,1));

    S_all = squeeze(nanstd(Cc_all,0,1));
    S_diff = squeeze(nanstd(Cc_diff,0,1));

    Tg_all = M_all./(S_all/sqrt(N_subj));
    Tg_diff = M_diff./(S_diff/sqrt(N_subj));

%     Tg_diff = squeeze(T_diff(c,:,:,:));
    
%      figure, hist(Tg_all(:),50);
%      figure, hist(S_all(:),50);
%     figure, hist(M_all(:),50);
%     figure, hist(S_bl(:),50);
%     figure, hist(Cc_all(:),50);
    
   
h = figure('visible',show_figures); h.Color='w'; pcolor(toi, foi, M_all./(S_all));  shading flat; colorbar;
    h = figure('visible',show_figures); h.Color='w'; pcolor(toi, foi, Tg_all);  shading flat; colorbar;
    title(sprintf('2nd-level t-value, Passing-Baseline [%s]', channels{c}));
    colormap(flipud(brewermap([],'RdBu')));
    saveas(h, sprintf('%s/timefreq_tval2_pass_bl_%s.png', figure_dir, channels{c}));
   
    h = figure('visible',show_figures); h.Color='w'; pcolor(toi, foi, Tg_diff); caxis([-6 6]);  shading flat; colorbar;
    title(sprintf('2nd-level t-value, Difficult-Easy [%s]', channels{c}));
    colormap(flipud(brewermap([],'RdBu')));
    saveas(h, sprintf('%s/timefreq_tval2_diff_easy_%s.png', figure_dir, channels{c}));
    
    df = N_subj - 1;
    Pg_all = 2 * (1 - tcdf(abs(Tg_all), df));
    Pg_diff = 2 * (1 - tcdf(abs(Tg_diff), df));
 
    p_thres = 0.05;
    p_str = '05';
    
    h = figure('visible',show_figures); h.Color='w'; pcolor(toi, foi, Tg_all.*single(Pg_all<p_thres)); caxis([-8 8]); shading flat; colorbar;
    title(sprintf('2nd-level t-value (p<%1.2f), Passing-Baseline [%s]', p_thres, channels{c}));
    colormap(flipud(brewermap([],'RdBu')));
    saveas(h, sprintf('%s/timefreq_tval2_p%s_all_%s.png', figure_dir, p_str, channels{c}));
    
    h = figure('visible',show_figures); h.Color='w'; pcolor(toi, foi, Tg_diff.*single(Pg_diff<p_thres)); caxis([-8 8]); shading flat; colorbar;
    title(sprintf('2nd-level t-value (p<%1.2f), Difficult-Easy [%s]', p_thres, channels{c}));
    colormap(flipud(brewermap([],'RdBu')));
    saveas(h, sprintf('%s/timefreq_pval2_p%s_diff_easy_%s.png', figure_dir, p_str, channels{c}));
    
    
end

%% Topo scalp plot

T_all = nan(length(channels_all), N_subj, length(foi), length(toi));
T_diff = nan(length(channels_all), N_subj, length(foi), length(toi));

for c = 1 : length(channels_all)
    Cc_all = squeeze(C_all(c,:,:,:));
    Cc_diff = squeeze(C_diff(c,:,:,:));
    
    M_all = squeeze(nanmean(Cc_all,1));
    M_diff = squeeze(nanmean(Cc_diff,1));

    S_all = squeeze(nanstd(Cc_all,0,1));
    S_diff = squeeze(nanstd(Cc_diff,0,1));

    Tg_all = M_all./(S_all/sqrt(N_subj));
    Tg_diff = M_diff./(S_diff/sqrt(N_subj)); 
    
    T_all(c,:,:,:) = Tg_all;
    T_diff(c,:,:,:) = Tg_diff;
    
end

