fft_all = squeeze(mean(data.fft.left_change.all.powspctrm,1));
idx = data.fft.left_change.all.freq<=50;
X = data.fft.left_change.all.freq(idx);

channels = [{'Fz'},{'FC5'},{'POz'},{'T7'}];
idx_c = [];
for i = 1 : length(channels)
    idx_c = [idx_c find(strcmp(eeg_channels,channels{i}))];
end

h = figure;
hh = plot(X, fft_all(idx_c,idx));

legend(channels);

fft_easy = squeeze(mean(data.fft.left_change.easy.powspctrm,1));
fft_diff = squeeze(mean(data.fft.left_change.difficult.powspctrm,1));

for i = 1 : length(channels)
    
    h = figure;
    Y = [fft_easy(idx_c(i),idx);fft_diff(idx_c(i),idx)]';
    
    hh = plot(X, Y);
    legend([{'Easy'},{'Difficult'}]);
    title(sprintf('FFT (Easy v Difficult) - %s', channels{i}));
    
    
end