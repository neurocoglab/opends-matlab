function [ data ] = run_ica( params, data )
% RUN_ICA Runs an ICA decomposition on the data and returns the result.
% Creates a new field in the "data" struct for the result, "data.eeg.ica"

cfg = params.eeg.ica.cfg;
cfg.channel = data.eeg.eeg_channels;

% Show console output for this one...
% data.eeg.ica = ft_componentanalysis(cfg, data.eeg.ft);
[~,data.eeg.ica] = evalc('ft_componentanalysis(cfg, data.eeg.ft);');

end

