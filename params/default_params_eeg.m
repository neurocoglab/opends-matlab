%% This script defines default EEG parameters for the opends processing 
%  pipeline

%% General stuff
params.eeg.bad_channel_file = '';
params.eeg.format = 'brainvision';
params.eeg.sub_dir = 'eeg';

%% Input stuff
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
params.eeg.ica.cfg = [];
params.eeg.ica.cfg.method = 'runica';
params.eeg.ica.plots.save = true;
params.eeg.ica.plots.cfg = [];
params.eeg.ica.plots.cfg.colormap = 'jet';
params.eeg.ica.plots.cfg.component = 1:20;
params.eeg.ica.plots.cfg.layout    = 'acticap-64ch-standard2.mat'; 
params.eeg.ica.plots.cfg.comment   = 'no';
params.eeg.ica.plots.window_size_topo = [1000 1000];
params.eeg.ica.plots.window_size_browser = [1500 900];

