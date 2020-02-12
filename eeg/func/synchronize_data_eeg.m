function [ data ] = synchronize_data_eeg( params, data )

% Get simulation time start for later synching
outdir = sprintf('%s/%s', params.io.output_dir, subject);
results_file = sprintf('%s/results_eye.mat',outdir);
flag_file = sprintf('%s/sim_logs.done', outdir);

if ~exist(flag_file, 'file')
   error('Simulation logs have not been preprocessed for subject "%s".\n', subject)
end

results_eye = load(results_file);
logid = 1;

t_simstart_sim = results_eye.data.sim.simstart.values{find(strcmp(results_eye.data.sim.events.hdr,'Millis'))}(1);
ids = results_eye.data.sim.events.values{find(strcmp(results_eye.data.sim.events.hdr,'LogId'))};
t_event_sim = results_eye.data.sim.events.values{find(strcmp(results_eye.data.sim.events.hdr,'Millis'))}(find(ids==logid));
t_delta = t_event_sim - t_simstart_sim; % Time between sim start and first event, in ms

data.eeg.t_delta_sim = t_delta;
data.eeg.t_start_sim = t_simstart_sim;

end

