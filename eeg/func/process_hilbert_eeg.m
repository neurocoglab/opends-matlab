function [ results ] = process_hilbert_eeg ( data, params, preproc, outdir )
%%%%%%%%%%%%%%
% Compute Hilbert transforms (amplitude envelopes) for specified frequency bands
% from the given EEG data
%


figdir = sprintf('%s/figures', outdir);

cfg_hb = data.eeg.cfg;
cfg_hb.continuous = 1;
cfg_hb.demean = 'no';

T = params.eeg.hilbert.bands;

N_bands = height(T);
results.eeg.hilbert.bands = T;

L = 0;
N_chan = length(data.eeg.ft.label);
for j = 1 : length(data.eeg.ft.trial)
    L = L + length(data.eeg.ft.time{j});
end

cfg_interp = [];
cfg_interp.method = 'linear';
cfg_interp.prewindow = 0.004;
cfg_interp.postwindow = 0.004;
cfg_interp.feedback = 'no';

[~,idata] = evalc('ft_interpolatenan_corrected(cfg_interp, data.eeg.ft);');
%idata = ft_interpolatenan_corrected(cfg_interp, data.eeg.ft);

for bb = 1 : N_bands

    freqband = [T.From(bb) T.To(bb)];
    fprintf(' Computing envelopes for %s: [%1.1f to %1.1f Hz]\n', ...
                T.Band{bb}, freqband(1), freqband(2));

    % make filter
    cfg_hb = [];
    cfg_hb.filtord = T.Filtord(bb);
    tol = 100;
    [b,a] = butter(cfg_hb.filtord,2*freqband/idata.fsample, 'bandpass');

    fmin = max(0.1,freqband(1)-10);
    fmax = max(20, freqband(2)+10);
    freqz(b,a,fmin:(fmax-fmin)/1000:fmax,idata.fsample);
    saveas(gcf,sprintf('%s/hilbert_filter_%s.png', figdir, T.Band{bb}));
    close(gcf);

    % apply filter to data
    N_t = size(idata.trial{1},2);
    datfilt = zeros(length(data.eeg.eeg_channels),N_t);
    envelopes = zeros(length(data.eeg.eeg_channels),N_t);
    for cc = 1 : length(data.eeg.eeg_channels)
        idx = find(strcmp(idata.label, data.eeg.eeg_channels{cc}));
        
        
        X = idata.trial{1}(idx,:);
        idx_nan = ~isnan(X);
        Y = filtfilt(b,a,X(idx_nan));
        datfilt(cc,idx_nan) = Y;
        envelopes(cc,idx_nan) = hilbert(Y);
    end

    results.eeg.hilbert.cfg(bb) = cfg_hb;
    results.eeg.hilbert.Fs = idata.fsample;
    results.eeg.hilbert.filtered(bb) = {datfilt};
    results.eeg.hilbert.channels = data.eeg.eeg_channels;
    results.eeg.hilbert.envelopes(bb) = {envelopes};

end

if params.eeg.hilbert.save_plots

    fprintf('Saving plots to %s.\n', figdir);
    
    et_results = load(sprintf('%s/results.mat', outdir)); 
    
    plot_params = preproc.params.plots.events;
    plot_params.plot_overtakes = true;
    plot_params.plot_saccades = true;
    plot_params.plot_lane_changes = true;
    plot_params.xlim=[0 2];

    plot_params.eeg.stdev = 6;
    plot_params.eeg.line_widths = [0.5 1.0 2];
    plot_params.eeg.alpha = [0.3, 0.6, 1.0];
    plot_params.eeg.height = 20;
    plot_params.eeg.interpolate = true;
    plot_params.eeg.smooth = 0;

    for bb = 1 : N_bands

       plot_params.to_file =  [{sprintf('%s/hilbert_ts_%s.fig', figdir, T.Band{bb})} ...
                               {sprintf('%s/hilbert_ts_%s.png', figdir, T.Band{bb})}];

       data_f = data.eeg.ft;
       data_f.trial={abs(results.eeg.hilbert.filtered{bb})};
       data_f.label = results.eeg.hilbert.channels;
       data_h = data.eeg.ft;
       data_h.trial={abs(results.eeg.hilbert.envelopes{bb})};
       data_h.label = results.eeg.hilbert.channels;
       
       plot_title = sprintf('Hilbert envelopes for %s [%1.1f to %1.1f Hz]', T.Band{bb}, T.From(bb), T.To(bb));
       plot_events2(et_results.results, plot_params, [{'Rounds'},{'LaneChanges'},{'Baseline'},{'Saccades'},{'Overtakes'}], ...
                    [{data.eeg.ft}, {data_f}, {data_h}], [{'AF7'},{'Fz'},{'FC5'},{'Cz'},{'T7'}], plot_title);
 
    end

end

end