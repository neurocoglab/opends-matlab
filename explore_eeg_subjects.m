%% Explore the EEG results of a single subject

subject = '0133';
% Get params
load processing_eeg_params
preproc = load('preproc_params_hd.mat');
proc = load('processing_params.mat');
ok = true;
outdir = sprintf('%s/%s/%s', preproc.params.root_dir, preproc.params.output_dir, subject);


%% ERPs - Difficult v. Easy

flag_file = sprintf('%s/eeg-erp.done', outdir);

if ~exist(flag_file,'file')
    fprintf('Subject %s has no data!', subject);
    ok = false;
end

if ok
   
    results_file = sprintf('%s/processing_results_eeg_erp.mat',outdir);
    load(results_file);
    
    % Average time series on scalp map
    
    figure;
    cfg = [];
    cfg.layout = 'acticap-64ch-standard2.mat';
    cfg.interactive = 'yes';
    cfg.showoutline = 'yes';
%     cfg.ylim=[-4 4];
    cfg.xlim=[-0.5 1.5];
    cfg.showlabels='yes';
    cfg.box='yes';

    h=ft_multiplotER(cfg, data_erp.left_change.easy, data_erp.left_change.difficult);
    
    clear data_erp;
    
end



%% FFT - Passing v. Baseline

if ok
    results_file = sprintf('%s/processing_results_eeg_fft.mat',outdir);
    load(results_file);
    
    norm_fft = false;
    show_channels = [{'Fz'},{'F8'},{'AF7'}];
    freq_range = [0 40];
    ampl_range = [0 6];
        
    figure;
    
    rows = 3;
    cols = ceil(length(show_channels) / rows);
    
    for i = 1 : length(show_channels)
        channel = show_channels{i};
        idx = find(strcmp(data_fft.epochs.baseline.label, channel));
        i1 = find(data_fft.epochs.baseline.freq > freq_range(1), 1);
        if isempty(i1), i1 = 1; end
        i2 = find(data_fft.epochs.baseline.freq > freq_range(2), 1);
        if isempty(i2), i2 = length(data_fft.epochs.baseline.freq); end
        idx_freq = i1 : i2;
        xb = data_fft.epochs.baseline.freq(idx_freq);
        Y = squeeze(data_fft.epochs.baseline.powspctrm(:,idx,idx_freq));
        Y = zscore(Y,0,2);
        fft_baseline = mean(Y,1);
        sterr_baseline = std(Y,0,1); % / ...
                         %size(data_fft.epochs.baseline.powspctrm,1);
        if norm_fft
            auc = trapz(xb,fft_baseline);
            fft_baseline = fft_baseline / auc;
        end
        subplot(rows, cols, i);
%         hh = plot(xb, fft_baseline);
        hh = plot_ci_filled(xb, fft_baseline, fft_baseline+sterr_baseline, fft_baseline-sterr_baseline, [0 0 1]);
        
%         hh.Color = 'b';
        hold on;
        
        idx = find(strcmp(data_fft.epochs.passing.label, channel));
        i1 = find(data_fft.epochs.passing.freq > freq_range(1), 1);
        if isempty(i1), i1 = 1; end
        i2 = find(data_fft.epochs.passing.freq > freq_range(2), 1);
        if isempty(i2), i2 = length(data_fft.epochs.passing.freq); end
        idx_freq = i1 : i2;
        xp = data_fft.epochs.passing.freq(idx_freq);
        Y = squeeze(data_fft.epochs.passing.powspctrm(:,idx,idx_freq));
        Y = zscore(Y,0,2);
        fft_passing = mean(Y,1);
        sterr_passing = std(Y,0,1);
        if norm_fft
            auc = trapz(xp,fft_passing);
            fft_passing = fft_passing / auc;
        end
        
        hh = plot_ci_filled(xp, fft_passing, fft_passing+sterr_passing, fft_passing-sterr_passing, [1 0 0]);
        hold on;
%         hh = plot(xp, fft_passing);
%         hh.Color = 'r';
        
        % Difference
        x1 = max(xb(1), xp(1));
        x2 = min(xb(end), xp(end));
        xx = x1:0.1:x2;
        yb = linterp(xb, fft_baseline, xx);
        yp = linterp(xp, fft_passing, xx);
        y_diff = yp - yb;
        
        hh = plot(xx, y_diff);
        hh.Color = 'k';
        hh.LineWidth = 2;
        hh.LineStyle = '--';
        
        title(sprintf('%s [%1.1f %1.1f]', channel, freq_range(1), freq_range(2)));
        legend([{'Baseline'},{'Passing'},{'P-D'}]);
        
        ylim(ampl_range);
        
    end
   
    
    
%     clear data_fft;
    
end

%% FFT - Difficult v. Easy

