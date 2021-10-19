function [ data ] = synchronize_eeg_sim( params, data, subject )

% Zero EEG time series on first simulation event trigger from the eye tracking
% data (note: this may not be a "simulation started" event)
%

outdir = sprintf('%s/%s', params.io.output_dir, subject);
results_file = sprintf('%s/results_preproc_eye.mat',outdir);
flag_file = sprintf('%s/flags/sim_logs_converted.done', outdir);

if ~exist(flag_file, 'file')
   error('Simulation logs have not been preprocessed for subject "%s".\n', subject)
end

T = load( results_file );
data.sim = T.data.sim;
clear T;

% Get EEG events, subtract time of zero event
if data.sim.zero_byte > 0
    idx0 = find(data.eeg.events.Trigger==data.sim.zero_byte, 1);
else
    error('\nNo zero serial byte defined for subject %s.', subject);
end
delta_t = data.eeg.events.Time(idx0);
data.eeg.events.Time = data.eeg.events.Time - delta_t;
data.eeg.ft.time = {data.eeg.ft.time{1} - delta_t};
data.eeg.t_simstart = delta_t;


end

