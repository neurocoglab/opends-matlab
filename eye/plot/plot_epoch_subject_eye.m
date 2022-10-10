function [ H ] = plot_epoch_subject_eye ( results, params, out2file )

    if nargin < 2
        out2file = false;
    end

    figdir = sprintf('%s/%s/figures', params.io.output_dir, results.subject);
    H = [];

    % Full session
    if out2file
        h = figure('visible','off');
    else
        h = figure;
    end
    H = [H h];
    set(h, 'Color', 'w');

    M = [results.eye.epochs.zscore.baseline.pupil;results.eye.epochs.zscore.nobaseline.pupil;results.eye.epochs.zscore.passing.pupil];
    grp = [zeros(length(results.eye.epochs.zscore.baseline.pupil),1);ones(length(results.eye.epochs.zscore.nobaseline.pupil),1); ...
                2*ones(length(results.eye.epochs.zscore.passing.pupil),1)];

    boxplot(M,grp,'outliersize',3,'symbol','');
    ylim(params.eye.epochs.plots.zlims);
    hh = title('Pupil diameter - Full session');
    set (hh,'FontSize',16);

    set(gca, 'XTickLabel', [{'Baseline'},{'Non-baseline'},{'Passing'}]);
    set (gca,'FontSize',14);

    hh = ylabel('Pupil diameter (z-score)');
    set (hh,'FontSize',14);

    if out2file
       saveas( h, sprintf('%s/epochs_all_eye.png', figdir) );
       close(h);
    end

    % By difficulty
    if params.sim.epochs.difficulty.apply
        
        if out2file
            h = figure('visible','off');
        else
            h = figure;
        end
        H = [H h];

        set(h, 'Color', 'w');

        M = [results.eye.epochs.zscore.passing_diff.pupil{1};results.eye.epochs.zscore.passing_diff.pupil{2}];
        grp = [ones(length(results.eye.epochs.passing_diff.pupil{1}),1)*results.eye.epochs.diff_levels(1); ...
               ones(length(results.eye.epochs.passing_diff.pupil{2}),1)*results.eye.epochs.diff_levels(2)];
        boxplot(M,grp,'outliersize',3,'symbol','');
        ylim(params.eye.epochs.plots.zlims);
        hh = title('Pupil diameter - Passing difficulty');
        set (hh,'FontSize',16);
        set(gca, 'XTickLabel', [{'Easy'},{'Difficult'}]);
        set (gca,'FontSize',14);
        hh = ylabel('Pupil diameter (z-score)');
        set (hh,'FontSize',14);

        if out2file
           saveas( h, sprintf('%s/epochs_difficulty_eye.png', figdir) );
           close(h);
        end
        
    end

    % Baselines over time
    if out2file
        h = figure('visible','off');
    else
        h = figure;
    end
    H = [H h];
    
    set(h, 'Color', 'w');

    N = length(results.eye.epochs.zscore.intervals.baseline.pupil);
    x = results.eye.epochs.intervals.baseline.times / 60000;
    y = zeros(N,1);
    err = zeros(N,1);

    for i = 1 : N
        pdi = results.eye.epochs.zscore.intervals.baseline.pupil{i};
        y(i) = mean(pdi);
        err(i) = std(pdi);  
    end

    colours = [[0 0 1];[1 0 0]];
    lighten = 0.8;

    hp = zeros(2, 1);
    hp(1) = plot_ci_filled( x, y, y+err, y-err, colours(1,:), lighten );
    %hh = errorbar(x,y,err,'b');
    hold on;

    N = length(results.eye.epochs.zscore.intervals.passing.pupil);
    x = results.eye.epochs.intervals.passing.times / 60000;
    y = zeros(N,1);
    err = zeros(N,1);

    for i = 1 : N
        pdi = results.eye.epochs.zscore.intervals.passing.pupil{i};
        y(i) = mean(pdi);
        err(i) = std(pdi);  
    end

    %errorbar(x,y,err,'r');
    hp(2) = plot_ci_filled( x, y, y+err, y-err, colours(2,:), lighten );

    hh = legend(hp,'Baseline','Passing');
    set (hh,'FontSize',14);

    % 
    hh = title('Pupil diameter - Epochs');
    set (hh,'FontSize',16);
    %        
    hh = ylabel('Pupil diameter (z-score)');
    set (hh,'FontSize',14);

    hh = xlabel('Time (min)');
    set (hh,'FontSize',14);
    
    ylim(params.eye.epochs.plots.zlims);

    if out2file
       saveas( h, sprintf('%s/epochs_time_eye.png', figdir) );
       close(h);
    end



end