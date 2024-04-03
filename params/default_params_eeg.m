%% This script defines default EEG parameters for the opends processing 
%  pipeline

%% General stuff
params.eeg.bad_channel_file = '';
params.eeg.sub_dir = 'eeg';

%% Input stuff
params.eeg.convert.format = 'brainvision';
params.eeg.convert.prefix = 'eeg-';
params.eeg.convert.suffix = '';
params.eeg.convert.downsample = [];
params.eeg.convert.restart_file = [];
params.eeg.convert.align_sim_triggers = false;
params.eeg.convert.start_byte = [];

params.eeg.layout = 'biosemi64.lay'; % 'acticap-64ch-standard2.mat'; 
params.eeg.electrodes = 'standard_1020.elc';

params.eeg.cfg = [];
params.eeg.cfg.dataset = '';
params.eeg.cfg.continuous = 'yes';

%% Preprocessing stuff

% Bandpass filter
params.eeg.bandpass.apply = true;
params.eeg.bandpass.cfg = [];
params.eeg.bandpass.cfg.lpfilter ='yes';
params.eeg.bandpass.cfg.hpfilter ='yes';
params.eeg.bandpass.cfg.bsfilter = 'no';
params.eeg.bandpass.cfg.bpfilter ='no';
params.eeg.bandpass.cfg.lpfreq = 60;
params.eeg.bandpass.cfg.hpfreq = 0.5;
params.eeg.bandpass.cfg.lpfiltord = 6;
params.eeg.bandpass.cfg.hpfiltord = 6;
params.eeg.bandpass.cfg.lpfilttype = 'but';
params.eeg.bandpass.cfg.hpfilttype = 'but';
params.eeg.bandpass.cfg.hpinstabilityfix = 'reduce';

% Notch
params.eeg.notch.apply = true;
params.eeg.notch.cfg = [];
params.eeg.notch.cfg.bsfilter = 'yes';
params.eeg.notch.cfg.bpfilter ='no';
params.eeg.notch.cfg.lpfilter ='no';
params.eeg.notch.cfg.hpfilter ='no';
params.eeg.notch.cfg.bsfreq = [48 62];
params.eeg.notch.cfg.bsfiltord = 4;
params.eeg.notch.cfg.bsfilttype = 'but';

% ICA
params.eeg.ica.apply = true;
params.eeg.ica.use_existing = false;
params.eeg.ica.cfg = [];
params.eeg.ica.cfg.method = 'runica';
params.eeg.ica.plots.save = true;
params.eeg.ica.plots.cfg = [];
params.eeg.ica.plots.cfg.colormap = 'jet';
params.eeg.ica.plots.cfg.component = 1:20;
params.eeg.ica.plots.cfg.layout    = params.eeg.layout;
params.eeg.ica.plots.cfg.comment   = 'no';
params.eeg.ica.plots.window_size_topo = [1000 1000];
params.eeg.ica.plots.window_size_browser = [1500 900];

% Artifact rejection
params.eeg.artifacts.plots.save = true;
params.eeg.artifacts.plots.channels = [{'vEOGover'},{'Fz'},{'FC5'},{'POz'},{'T7'}];
params.eeg.artifacts.plots.stdev = 6;
params.eeg.artifacts.eye.apply = true;

params.eeg.artifacts.zscore.apply = true;
cfg = [];
cfg.artfctdef.zvalue.cutoff = 20;
cfg.artfctdef.zvalue.trlpadding = 0;
cfg.artfctdef.zvalue.fltpadding = 0;
cfg.artfctdef.zvalue.artpadding = 0.05;
params.eeg.artifacts.zscore.cfg = cfg;

params.eeg.artifacts.minaccepttim = 0.7;

params.eeg.artifacts.remove_channels = {};

% Interpolation
params.eeg.artifacts.interpolate.apply = true;
cfg = [];
cfg.method = 'linear';
cfg.prewindow = 0.004;
cfg.postwindow = 0.004;
params.eeg.artifacts.interpolate.cfg = cfg;



%% Processing stuff

% Trials
params.eeg.trials.events = [{'LeftChange'}, ...
                            {'RightChange'}, ...
                            {'Brake'}];

params.eeg.trials.windows_file = 'trial_windows_eeg.csv';
                        

% Hilbert analysis
params.eeg.hilbert.apply = true;
params.eeg.hilbert.bands_file = 'hilbert_bands_eeg.csv';
params.eeg.hilbert.downsample = 250;
params.eeg.hilbert.save_filtered = false;

params.eeg.hilbert.outlier_zscore = 3;

params.eeg.hilbert.plots.save = true;
params.eeg.hilbert.plot.show_plots = true;
params.eeg.hilbert.plot.bands = 'all';
params.eeg.hilbert.plot.channels = 'all';
params.eeg.hilbert.plot.scale_pupil = 2;
params.eeg.hilbert.plot.smooth_pupil = 70;
params.eeg.hilbert.plot.scale_envelope = 1;
params.eeg.hilbert.plot.scale_general = 0.2;
params.eeg.hilbert.plot.xlims = [0 4];
params.eeg.hilbert.plot.topo.baseline_overtake.clim = [-3 3];
params.eeg.hilbert.plot.topo.baseline_overtake.colormap = 'RdBu';

params.eeg.hilbert.analysis.p_crit = 0.01;
params.eeg.hilbert.analysis.p_crit_fdr = 0.05;

% ERP analysis
params.eeg.erp.apply = true;


% Time/frequency analysis
params.eeg.timefreq.apply = true;
params.eeg.timefreq.buffer = 5;
params.eeg.timefreq.event_types_file = 'timefreq_events_eeg.csv';
params.eeg.timefreq.dt = 0.1;
params.eeg.timefreq.foi = [2:1:8,10:2:18,22:4:30];
params.eeg.timefreq.wavelet = 'cmor4-1';
params.eeg.timefreq.regen_all = false;
