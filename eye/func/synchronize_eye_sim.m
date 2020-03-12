function [ data ] = synchronize_eye_sim( params, data )
%SYNCHRONIZE_EYE_SIM Synchronize eye and simulation time series to the
%first event trigger

sim_dir = sprintf('%s/%s/%s', params.io.output_dir, data.subject, params.sim.sub_dir);
input_file = sprintf('%s/events-All.csv', sim_dir);
data.sim.events.values = import_log_sim(input_file, params.sim.log.events_format);

log_ids = data.eye.log.messages.LogId;
log_times = data.eye.log.messages.Time;

% Zero eye time on first trigger
data.eye.t_start_sim = double(log_times(1));
data.eye.t = data.eye.t - data.eye.t_start_sim;

% Zero sim time on first trigger (note: this may not be the simulation
% start trigger as this isn't inserted in older versions)
data.sim.zero_id = log_ids(1);
idx_sim = find(data.sim.events.values.LogId==log_ids(1), 1);

if isempty(idx_sim)
   error('\tNo initial event (log id = %d) found in simulation log!', log_ids(1)); 
end

data.sim.zero_byte = data.sim.events.values.AdjSerialByte(idx_sim);
data.sim.t_start = data.sim.events.values.Millis(idx_sim);


end

