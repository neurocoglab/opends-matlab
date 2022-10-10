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
Fs = idata.fsample;

for bb = 1 : N_bands

    freqband = [T.From(bb) T.To(bb)];
    fprintf('\tComputing envelopes for %s: [%1.1f to %1.1f Hz]', ...
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
    
    Fs2 = Fs;
    %N_t = size(idata.trial{1},2);
    if params.eeg.hilbert.save_filtered
        datfilt = cell(length(eeg_channels),1);
    end
    envelopes = cell(length(eeg_channels),1);
    times = cell(length(eeg_channels),1);
    for cc = 1 : length(eeg_channels)
        idx = strcmp(idata.label, eeg_channels{cc});
        
        X = idata.trial{1}(idx,:);
        idx_nan = ~isnan(X);
        Y = filtfilt(b,a,X(idx_nan));
        time_cc = idata.time{1}(idx_nan);
        
        datfilt_cc = Y;
        H = hilbert(Y);
        
         % Downsample?
        if params.eeg.hilbert.downsample > 0 && params.eeg.hilbert.downsample < Fs
            step = round(idata.fsample/params.eeg.hilbert.downsample);
            ts_env = timeseries(H, time_cc);
            idx_ds = 1:step:length(ts_env.Time);
            ts_env = resample(ts_env, ts_env.Time(idx_ds));
            H = ts_env.Data;
            if params.eeg.hilbert.save_filtered
                ts_filt = timeseries(datfilt_cc, time_cc);
                ts_filt = resample(ts_filt, ts_filt.Time(idx_ds));
                datfilt_cc = ts_filt.Data;
            end
            Fs2 = 1/(ts_env.Time(2)- ts_env.Time(1)); % Time is in seconds
            time_cc = single(ts_env.Time);
        end
        
        envelopes(cc) = {squeeze(single(H))};
        times(cc) = {time_cc};
        if params.eeg.hilbert.save_filtered
            datfilt(cc) = {single(datfilt_cc)};
        end
    end
    
    if params.general.debug
        fprintf(' [downsampled from %1.2f to %1.2f Hz]', Fs, Fs2);
    end
            
    results.eeg.hilbert.cfg(bb) = cfg_hb;
    results.eeg.hilbert.Fs = Fs2;
    if params.eeg.hilbert.save_filtered
        results.eeg.hilbert.filtered(bb) = {datfilt};
    end
    results.eeg.hilbert.channels = eeg_channels;
    results.eeg.hilbert.envelopes(bb) = {envelopes};
    results.eeg.hilbert.time(bb) = {times};
    
    fprintf('\n');

end

end