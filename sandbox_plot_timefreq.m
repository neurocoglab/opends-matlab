%% Plot time/freq with axes

close all;

fig_dir = '/Users/lpzatr/OneDrive/OneDrive UoN/OneDrive - The University of Nottingham/synched/projects/driving/figures';

trial_diff = [{'Easy'},{'Difficult'}];
which_y = [{'easy'},{'diff'}];
chans = [{'Fz'},{'F8'},{'AF7'}];

X = [fftdata.all{1}];

idx_f = find(X.freq < 20 & X.freq > 2);

for d = 1 : 2
    
    Y = eval(sprintf('Y_%s;', which_y{d}));
    X.powspctrm = squeeze(nanmean(Y,1));

    for ch = 1 : length(chans)

        idx_c = find(strcmp(X.label, chans{ch}));
        P = squeeze(X.powspctrm(idx_c,idx_f,:));

        h = figure;
        h.Color = 'w';
        hh = imagesc(flip(P,1));

        hh = title(sprintf('%s Trials: %s', trial_diff{d}, chans{ch}));
        hh.FontSize = 26;

        colormap('jet');
        caxis([0 6000]);

        ax = gca;

        F = flip(X.freq(idx_f));
        idx_t = 1:floor(length(F) / 10):length(F);
        labels = [];
        for i = 1 : length(idx_t)
            labels = [labels {sprintf('%1.1f', F(idx_t(i)))}];
        end
        ax.YAxis.TickValues = idx_t;
        ax.YAxis.TickLabels = labels;

        hh = ylabel('Frequency (Hz)');
        hh.FontSize = 16;

        T = X.time;
        tstep = 0.2;
        step = round((T(end)-T(1))/tstep);
        idx_0 = find(T==0);
        idx_1 = flip(idx_0:-step:1);
        idx_2 = idx_0+step:step:length(T);
        idx_t = [idx_1 idx_2];
        labels = [];
        for i = 1 : length(idx_t)
            labels = [labels {sprintf('%1.1f', T(idx_t(i)))}];
        end
        ax.XAxis.TickValues = idx_t;
        ax.XAxis.TickLabels = labels;

        hh = xlabel('Time from left-change (s)');
        hh.FontSize = 16;

        saveas(h, sprintf('%s/timefreq_%s_%s.png', fig_dir, which_y{d}, chans{ch}));
        
    end


end