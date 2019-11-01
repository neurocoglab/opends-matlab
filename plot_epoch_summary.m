function [ h ] = plot_epoch_summary ( summary, subject_data, params, out2file )

if nargin < 4
   out2file = false; 
end

% Full session (z-score)
if out2file
    h = figure('visible','off');
    outdir = sprintf('%s/figures', params.output_dir);
    if ~exist(outdir, 'dir')
       mkdir(outdir); 
    end
else
    h = figure;
end

set(h,'Color','w');

M = summary.stats.passing_baseline.data.PD;
M = [summary.zscore.baseline.pupil;summary.zscore.passing.pupil];
grp = [zeros(length(summary.zscore.baseline.pupil),1);ones(length(summary.zscore.passing.pupil),1)];

boxplot(M,grp,'outliersize',3,'symbol','');

ylim([-5 5]);

hh = title('Pupil diameter - All subjects');
set (hh,'FontSize',16);
       
set(gca, 'XTickLabel', [{'Baseline'},{'Passing'}]);
set (gca,'FontSize',14);

hh = ylabel('Pupil diameter (z-score)');
set (hh,'FontSize',14);

resize_window(h, [1000,600], [500,500]);

if out2file
%     saveas(h, sprintf('%s/lum_scatter.fig', outdir));
    saveas(h, sprintf('%s/boxplots_pd_zscore_pass_baseline.png', outdir));
    close(h);
end


% By difficulty
if out2file
    h = figure('visible','off');
else
    h = figure;
end

set(h,'Color','w');
M = [];
grp = [];

for i = 1 : length(summary.diff_levels)
   x = summary.zscore.passing_diff.pupil{i};
   M = [M;x];
   x = i*ones(length(x),1);
   grp = [grp;x];
end

boxplot(M,grp,'outliersize',3,'symbol','');
hh = title('Pupil diameter - Passing difficulty');
set (hh,'FontSize',16);
% set(gca, 'XTickLabel', summary.diff_levels);
set(gca, 'XTickLabel', [{'Easy'},{'Difficult'}]);
set (gca,'FontSize',14);
hh = ylabel('Pupil diameter (z-score)');
set (hh,'FontSize',14);

ylim([-3 4]);

resize_window(h, [1000,600], [500,500]);

if out2file
    saveas(h, sprintf('%s/boxplots_pdz_pass_difficulty.png', outdir));
    close(h);
end

% By outcome
if out2file
    h = figure('visible','off');
else
    h = figure;
end

set(h,'Color','w');

M = [summary.zscore.passing_outcome.positive.pupil;summary.zscore.passing_outcome.negative.pupil];
grp = [zeros(length(summary.zscore.passing_outcome.positive.pupil),1);ones(length(summary.zscore.passing_outcome.negative.pupil),1)];

boxplot(M,grp,'outliersize',3,'symbol','');
hh = title('Pupil diameter - Outcome valence');
set (hh,'FontSize',16);
% set(gca, 'XTickLabel', summary.diff_levels);
set(gca, 'XTickLabel', [{'Positive'},{'Negative'}]);
set (gca,'FontSize',14);
hh = ylabel('Pupil diameter (z-score)');
set (hh,'FontSize',14);

ylim([-3 4]);

resize_window(h, [1000,600], [500,500]);

if out2file
    saveas(h, sprintf('%s/boxplots_pdz_pass_outcome.png', outdir));
    close(h);
end

% Saccade rate

% Full session (raw)
if out2file
    h = figure('visible','off');
else
    h = figure;
end

set(h,'Color','w');

M = [summary.baseline.saccade_rate;summary.passing.saccade_rate];
grp = [zeros(length(summary.baseline.saccade_rate),1);ones(length(summary.passing.saccade_rate),1)];

boxplot(M,grp,'outliersize',3,'symbol','');

ylim([-0.5 3]);

hh = title('Saccade rate - All subjects');
set (hh,'FontSize',16);
       
set(gca, 'XTickLabel', [{'Baseline'},{'Passing'}]);
set (gca,'FontSize',14);

hh = ylabel('Saccade rate (a.u.)');
set (hh,'FontSize',14);

resize_window(h, [1000,600], [500,500]);

if out2file
    saveas(h, sprintf('%s/boxplots_saccade_pass_baseline.png', outdir));
    close(h);
end

if out2file || params.epochs.plot_scatter

    % Subject score versus pupil dilation
    sd_subs = subject_data.UID;
    y_baseline = zeros(length(summary.subjects),1);
    y_passing = zeros(length(summary.subjects),1);
    diff_bp = zeros(length(summary.subjects),1);
    x = zeros(length(summary.subjects),1);

    clr = [0.1 1 1];

    X = subject_data.Game_score;

    for i = 1 : length(summary.subjects)

        subject = summary.subjects{i};
        idx = find(strcmp(sd_subs, subject));

        y_baseline(i) = mean(summary.zscore.baseline.subjects.pupil{i});
        y_passing(i) = mean(summary.zscore.passing.subjects.pupil{i});
        diff_bp(i) = y_passing(i) - y_baseline(i);

        x(i) = X(idx(1));

    end

    mi=min(x)-200; mx = max(x);
    
    if out2file
        h = figure('visible','off');
    else
        h = figure;
    end
    
    set(h,'Color','w');
    offset = 1:length(summary.subjects);

    for i = 1 : length(summary.subjects)
        cc = 1 - clr * (x(i)-mi)/(mx-mi);
        plot([1 2],[y_baseline(i)+offset y_passing(i)+offset], 'o-', 'Color', cc, 'LineWidth', 1.5);
        hold on;
    end

    set(gca,'XTickLabel',[{''},{'Baseline'},{''},{'Passing'},{''}]);
    title('Baseline-Passing Differences');
    xlim([0.5 2.5]);
    ylabel('Pupil Diameter (Z score)');
    
    resize_window(h, [1000,600], [500,500]);

    if out2file
        saveas(h, sprintf('%s/lines_subjects_pass_baseline.png', outdir));
        close(h);
    end

    % Scatterplots (raw)
%     if out2file
%         h = figure('visible','off');
%     else
%         h = figure;
%     end
%     
%     scatter(x,y_baseline,'b');
%     hh=lsline(gca);
%     set(hh,'Color','b');
%     set(hh,'LineWidth',1.5);
%     set(h,'Color','w');
%     title('Baseline PD / Score');
%     xlabel('Simulation Score');
%     ylabel('Mean pupil diameter (mm)');
%     [r,p] = corr(x,y_baseline);
%     fprintf('Baseline correlation: %1.4f (p=%1.4f)\n', r, p);
%     
%     if out2file
%         saveas(h, sprintf('%s/lines_subjects_pass_baseline.png', outdir));
%         close(h);
%     end
% 
%     h=figure; 
%     
%     scatter(x,y_passing,'r');
%     hh=lsline(gca);
%     set(hh,'Color','r');
%     set(hh,'LineWidth',1.5);
%     set(h,'Color','w');
%     title('Baseline PD / Score');
%     xlabel('Simulation Score');
%     ylabel('Mean pupil diameter (mm)');
% 
%     title('Passing PD / Score');
%     [r,p] = corr(x,y_passing);
%     fprintf('Passing correlation: %1.4f (p=%1.4f)\n', r, p);
%     if params.epochs.save_plots
%        print(sprintf('%s/%spassing_score_scatter_raw.png',params.output_dir,params.epochs.plot_prefix),'-dpng');
%     end
% 
%     h=figure; scatter(x,diff_bp,'k');
%     hh=lsline(gca);
%     set(hh,'Color','k');
%     set(hh,'LineWidth',1.5);
%     set(h,'Color','w');
%     title('Baseline PD / Score');
%     xlabel('Simulation Score');
%     ylabel('Mean pupil diameter (mm)');
%     title('Baseline - PD / Score');
%     [r,p] = corr(x,diff_bp);
%     fprintf('P-B correlation: %1.4f (p=%1.4f)\n', r, p);
%     if params.epochs.save_plots
%        print(sprintf('%s/%sbp_score_scatter_raw.png',params.output_dir,params.epochs.plot_prefix),'-dpng');
%     end

    % Scatterplots (z-score)
    y_baseline = zeros(length(summary.subjects),1);
    y_passing = zeros(length(summary.subjects),1);
    diff_bp = zeros(length(summary.subjects),1);

    for i = 1 : length(summary.subjects)
        y_baseline(i) = mean(summary.zscore.baseline.subjects.pupil{i});
        y_passing(i) = mean(summary.zscore.passing.subjects.pupil{i});
        diff_bp(i) = y_passing(i) - y_baseline(i);

    end

    if out2file
        h = figure('visible','off');
    else
        h = figure;
    end
    
    set(h,'Color','w');
    
    scatter(x,y_baseline,'b');
    hh=lsline(gca);
    set(hh,'Color','b');
    set(hh,'LineWidth',1.5);
    hold on;    
    scatter(x,y_passing,'r');
    hh=lsline(gca);
    set(hh,'Color','r');
    set(hh,'LineWidth',1.5);
    
    title('Pupil diameter X Score');
    xlabel('Simulation Score');
    ylabel('Mean pupil diameter (z-score)');
    
    [r,p] = corr(x,y_baseline);
    fprintf('Baseline correlation: %1.4f (p=%1.4f)\n', r, p);
    [r,p] = corr(x,y_passing);
    fprintf('Passing correlation: %1.4f (p=%1.4f)\n', r, p);
    
    resize_window(h, [1000,600], [500,500]);
    
    if out2file
        saveas(h, sprintf('%s/scatter_pdz_X_score.png', outdir));
        close(h);
    end

    if out2file
        h = figure('visible','off');
    else
        h = figure;
    end 
    
    scatter(x,diff_bp,'k');
    hh=lsline(gca);
    set(hh,'Color','k');
    set(hh,'LineWidth',1.5);
    set(h,'Color','w');
    title('Baseline PD / Score');
    xlabel('Simulation Score');
    ylabel('Mean pupil diameter (z-score)');
    title('Baseline - PD / Score');
    [r,p] = corr(x,diff_bp);
    fprintf('P-B correlation: %1.4f (p=%1.4f)\n', r, p);
    
    if out2file
        saveas(h, sprintf('%s/scatter_pdz_bl-pd_X_score.png', outdir));
        close(h);
    end

end

% Plot cycles
N_cycles = length(summary.cycles.baseline.pupil);

% Baseline - Boxplots


% Passing - Boxplots


% Baseline & Passing - Line plot

pd_bl_mean = zeros(N_cycles,1);
pd_bl_ci = zeros(N_cycles,2);
pd_pass_mean = zeros(N_cycles,1);
pd_pass_ci = zeros(N_cycles,2);

N_subs = length(summary.subjects);

pd_bl = zeros(N_cycles, N_subs);
pd_pass = zeros(N_cycles, N_subs);

for j = 1 : N_cycles
    x_bl = zeros(N_subs,1);
    x_pass = zeros(N_subs,1);
    for i = 1 : N_subs
       if length(summary.zscore.cycles.baseline.subjects.pupil{i}) < N_cycles
           %ffs
           x_bl(i) = nan;
           x_pass(i) = nan;
       else
           x_bl(i) = mean(summary.zscore.cycles.baseline.subjects.pupil{i}{j});
           x_pass(i) = mean(summary.zscore.cycles.passing.subjects.pupil{i}{j});
           
       end
    end
    
    pd_bl(j,:) = x_bl;
    pd_pass(j,:) = x_pass;
    
    pd_bl_mean(j) = nanmean(x_bl);
    pd_pass_mean(j) = nanmean(x_pass);
    
    pdist = fitdist(x_bl(~isnan(x_bl)),'Normal');
    ci = paramci(pdist);
    pd_bl_ci(j,:) = ci(:,1);
    pdist = fitdist(x_pass(~isnan(x_pass)),'Normal');
    ci = paramci(pdist);
    pd_pass_ci(j,:) = ci(:,1);
    
end

% Regress
% df1 = 1; df2 = N_subs-1;
% X=[ones(N_subs*N_cycles,1),repmat(1:N_cycles,1,N_subs)'];
% 
% y=reshape(pd_bl,N_subs*N_cycles,1);
% [b,bint,r,rint,stats] = regress(y,X);
% fprintf('Regression - Baseline PD ~ Round: b=%1.3f, r2=%1.3f, F(%d,%d)=%1.3f; p=%1.3f\n', ...
%                                         b(2), stats(1), df1, df2, stats(2), stats(3));
% 
% y=reshape(pd_pass,N_subs*N_cycles,1);
% [b,bint,r,rint,stats] = regress(y,X);
% fprintf('Regression - Passing PD ~ Round: b=%1.3f, r2=%1.3f, F(%d,%d)=%1.3f; p=%1.3f\n', ...
%                                         b(2), stats(1), df1, df2, stats(2), stats(3));
%                                     
% % T-tests
% tvals = zeros(N_cycles,1);
% pvals = tvals;
% statstr = cell(N_cycles,1);
% for j = 1 : N_cycles
%     [h,p,ci,stats] = ttest(pd_bl(j,:),pd_pass(j,:));
%     pvals(j) = p;
%     tvals(j) = stats.tstat;
%     ast='ns';
%     if p<0.001
%         ast = '**';
%     elseif p < 0.01   
%         ast = '*';
%     end
%     statstr(j) = {sprintf('%1.2f(%s) ', stats.tstat, ast)};
% end
% 
% fprintf(['T-values (PD baseline versus passing): ' [statstr{:}] '\n']);
                                    
if out2file
    h = figure('visible','off');
else
    h = figure;
end

set(h, 'Color', 'w');

x = [1:N_cycles; 1:N_cycles]';
y = [pd_bl_mean,pd_pass_mean];
lower = [pd_bl_ci(:,1),pd_pass_ci(:,1)]; % / N_subs;
upper = [pd_bl_ci(:,2),pd_pass_ci(:,2)];
colours = [[0.8 0 0];[0 0.8 0.2]];
plot_ci_filled(x, y, upper, lower, colours, 0.8, '-o');

ylim([-1.5 2]);
legend([{'Baseline'},{'Passing'}]);

xlabel('Simulation Round');
ylabel('Mean pupil diameter (z-score)');

set (gca,'FontSize',14);
hh = title('Pupil diameter - Over Time');
set (hh,'FontSize',16);

resize_window(h, [1000,600], [500,500]);

if out2file
    saveas(h, sprintf('%s/lines_pdz_rounds.png', outdir));
    close(h);
end

% Boxplot
N_min = Inf;
for i = 1 : N_cycles
   N_min = min(N_min, length(summary.zscore.cycles.baseline.pupil{i}));
   N_min = min(N_min, length(summary.zscore.cycles.passing.pupil{i}));
end

if out2file
    h = figure('visible','off');
else
    h = figure;
end

set(h,'Color','w');

data = zeros(N_subs, N_cycles, 2);
data(:,:,1) = pd_bl';
data(:,:,2) = pd_pass';

data(isnan(data))=mean(data(:));

iosr.statistics.boxPlot(data,'theme','colorall','showViolin',false,'showOutliers',false, ...
                             'groupLabels',[{'Baseline'},{'Passing'}],'showLegend',true);

xlabel('Simulation Round');
ylabel('Mean pupil diameter (z-score)');

set (gca,'FontSize',14);
hh = title('Pupil diameter - Over Time');
set (hh,'FontSize',16);

resize_window(h, [1000,600], [500,500]);

if out2file
    saveas(h, sprintf('%s/bxplots_pdz_rounds.png', outdir));
    close(h);
end



end