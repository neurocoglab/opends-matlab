if params.eeg.timefreq.plots.plot_rounds
    h=figure('visible','on');
    h.Color = 'w';
    for r = 1 : 8
        tfs = results.eeg.timefreq.left_change.rounds{r}.timelock; 
        subplot(4,2,r),pcolor(tfs.time,tfs.freq,squeeze(tfs.avg(38,:,:)));caxis([0 3]);shading flat;
        title(sprintf('Round %d', r));
    end
    suptitle('Left Change by Round - Fz');
    h.Position=[0 0 1000 800];
%     saveas(h,
end

h=figure('visible','on');
h.Color = 'w';

tfs = results.eeg.timefreq.left_change.easy.timelock;
subplot(1,2,1),pcolor(tfs.time,tfs.freq,squeeze(tfs.avg(38,:,:)));caxis([0 3]);shading flat;
title('Easy');
colorbar;
tfs = results.eeg.timefreq.left_change.difficult.timelock;
subplot(1,2,2),pcolor(tfs.time,tfs.freq,squeeze(tfs.avg(38,:,:)));caxis([0 3]);shading flat;
title('Difficult');
colorbar;

suptitle('Left Change by Difficulty - Fz');
h.Position=[220 220 1200 500];

h=figure('visible','on');
h.Color = 'w';

tfs = results.eeg.timefreq.left_change.negative.timelock;
subplot(1,2,1),pcolor(tfs.time,tfs.freq,squeeze(tfs.avg(38,:,:)));caxis([0 3]);shading flat;
title('Negative');
colorbar;
tfs = results.eeg.timefreq.left_change.positive.timelock;
subplot(1,2,2),pcolor(tfs.time,tfs.freq,squeeze(tfs.avg(38,:,:)));caxis([0 3]);shading flat;
title('Positive');
colorbar;

suptitle('Left Change by Outcome - Fz');
h.Position=[230 230 1200 500];