function [ results ] = process_hilbert_eeg ( params, data, results )
%%%%%%%%%%%%%%
% Compute Hilbert transforms (amplitude envelopes) for specified frequency bands
% from the given EEG data
%

if ~exist('results', 'var') 
    results = [];
end

results.subject = data.subject;

outdir = sprintf( '%s/%s', params.io.output_dir, data.subject );
figdir = sprintf( '%s/figures', outdir );

cfg_hb = data.eeg.cfg;
cfg_hb.continuous = 1;
cfg_hb.demean = 'no';

% Load Hilbert bands and filter orders from CSV file
fn = params.eeg.hilbert.bands_file;
T = readtable( fn );

N_bands = height(T);
results.eeg.hilbert.bands = T;

idata = data.eeg.ft;

for bb = 1 : N_bands

    freqband = [T.From(bb) T.To(bb)];
    fprintf('\tComputing envelopes for %s: [%1.1f to %1.1f Hz]\n', ...
                T.Band{bb}, freqband(1), freqband(2));

    % make filter
    cfg_hb = [];
    cfg_hb.filtord = T.Filtord(bb);
    [b,a] = butter(cfg_hb.filtord,2*freqband/idata.fsample, 'bandpass');

    if params.eeg.hilbert.plots.save
        fmin = max(0.1,freqband(1)-10);
        fmax = max(20, freqband(2)+10);
        freqz(b,a,fmin:(fmax-fmin)/1000:fmax,idata.fsample);
        saveas(gcf,sprintf('%s/eeg_hilbert_filter_%s.png', figdir, T.Band{bb}));
        close(gcf);
    end

    % apply filter to data
    eeg_channels = data.eeg.eeg_channels;
    for k = 1 : length(data.eeg.bad_channels)
        eeg_channels(strcmp(eeg_channels, data.eeg.bad_channels{k}))=[];
    end
    
    N_t = size(idata.trial{1},2);
    datfilt = zeros(length(eeg_channels),N_t);
    envelopes = zeros(length(eeg_channels),N_t);
    for cc = 1 : length(eeg_channels)
        idx = find(strcmp(idata.label, eeg_channels{cc}));
        
        X = idata.trial{1}(idx,:);
        idx_nan = ~isnan(X);
        Y = filtfilt(b,a,X(idx_nan));
        datfilt(cc,idx_nan) = Y;
        envelopes(cc,idx_nan) = hilbert(Y);
    end

    results.eeg.hilbert.cfg(bb) = cfg_hb;
    results.eeg.hilbert.Fs = idata.fsample;
    results.eeg.hilbert.filtered(bb) = {single(datfilt)};
    results.eeg.hilbert.channels = eeg_channels;
    results.eeg.hilbert.envelopes(bb) = {single(envelopes)};
    results.eeg.hilbert.time = data.eeg.ft.time{1};

end

end