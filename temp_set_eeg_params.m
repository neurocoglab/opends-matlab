clear params

params.eeg.data_dir = '/Volumes/AndrewElements/data/driving/eeg';

params.eeg.cfg = [];
params.eeg.cfg.headerfile = [];
params.eeg.cfg.dataset = [];
params.eeg.cfg.continuous = true;
params.eeg.cfg.blocksize = 60;

params.eeg.cfg.hpfilter = 'no';
params.eeg.cfg.hpfreq = 1;
params.eeg.cfg.lpfilter = 'no';
params.eeg.cfg.lpfreq = 100;

