%%%%
% Computes the time offset require to convert the eye tracker time variable
% from time relative to its first timepoint, to time relative to the
% simulation start; this allows time series in tracker and EEG to be
% aligned. Adding the output (t_delta) will align the time series.
%
% Returns:  t_delta (ms)

function [t_delta, t_simstart_sim] = get_simstart_offset( params, subject )

    outdir = sprintf('%s/%s/%s', params.root_dir, params.output_dir, subject);
    results_file = sprintf('%s/results.mat',outdir);
    flag_file = sprintf('%s/sim_logs.done', outdir);

    if ~exist(flag_file, 'file')
       error('ET logs have not been preprocessed for subject "%s".\n', subject)
    end

    load(results_file);

    % Get first event time [tracker time AND log time]; get sim start time [log time]
    % Compute sim start time as tracker time; subtract from time series

    logid = 1;

    t_simstart_sim = data.sim.simstart.values{find(strcmp(data.sim.events.hdr,'Millis'))}(1);
    ids = data.sim.events.values{find(strcmp(data.sim.events.hdr,'LogId'))};
    t_event_sim = data.sim.events.values{find(strcmp(data.sim.events.hdr,'Millis'))}(find(ids==logid));
    t_delta = t_event_sim - t_simstart_sim; % Time between sim start and first event, in ms
    ids = data.eye.log.messages{find(strcmp(data.eye.log.hdr,'LogId'))};
    t_event_tracker =  data.eye.log.messages{find(strcmp(data.eye.log.hdr,'Time'))}(find(ids==logid));

    t_event_tracker = int64(t_event_tracker / 1000); % tracker unit is microseconds, convert to ms
    t_simstart_tracker = t_event_tracker - t_delta; % simulation start in tracker time

    t_delta = single(t_simstart_tracker - data.eye.t_start); % Time is currently relative to first timepoint,
                                                             % make it relative to sim start 


end
