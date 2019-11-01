function [ data, results ] = process_rounds( data, results, params, subject )
% PROCESS_ROUNDS - Process the simulation rounds using sim events,
%                  mapping simulation time to eye tracker time. Also
%                  maps baseline distances to eye tracker time series.

% Load blink-processed data
% infile = sprintf('%s/saccades.mat',output_dir);
% load(infile);

% Does this subject have a lane distance specified? (Fix for bug in early
% version)
lane_dist = -1;
n_rounds = -1;

if ~isempty(params.rounds.fixdist)
    subjs = params.rounds.fixdist{1};
    idx = find(strcmp(subjs, subject));
    if ~isempty(idx)
       dists = params.rounds.fixdist{2};
       lane_dist = dists(idx);
       rounds = params.rounds.fixdist{3};
       n_rounds = rounds(idx);
    end
end

output_dir = sprintf('%s/%s', params.root_dir, params.output_dir);

% Load messages from eye tracker
input_file = sprintf('%s/%s/%s_messages.csv', output_dir, subject, subject);

if exist(input_file, 'file')
   [data.eye.log.messages, data.eye.log.hdr] = import_messages(input_file);
else
   % There are probably multiple parts; need to glue these together
   input_file = sprintf('%s/%s/%s_part1_messages.csv', output_dir, subject, subject);
   if ~exist(input_file, 'file')
       error('No message log found for subject %s', subject);
   end
   
   messages = [];
   k = 1;
   while exist(input_file, 'file')
       [messages_i, hdr] = import_messages(input_file);
       if k == 1
           messages = messages_i;
       else
           for j = 1 : 4
               messages(j) = {[messages{j};messages_i{j}]};
           end
       end
       
       
       k = k + 1;
       input_file = sprintf('%s/%s/%s_part%d_messages.csv', output_dir, subject, subject, k);
       
   end
   
   data.eye.log.messages = messages;
   data.eye.log.hdr = hdr;
end

% Load TriggerActuated | LaneChange | IncreaseSpeed events from simulation log
input_file = sprintf('%s/%s/events-All.csv', output_dir, subject);
[data.sim.events.values, data.sim.events.hdr] = import_log(input_file, params.simlog.events_format);
input_file = sprintf('%s/%s/events-SimulatorStarted.csv', output_dir, subject);
[data.sim.simstart.values, data.sim.simstart.hdr] = import_log(input_file, params.simlog.start_format);
input_file = sprintf('%s/%s/events-TriggerActuated.csv', output_dir, subject);
[data.sim.triggers.values, data.sim.triggers.hdr] = import_log(input_file, params.simlog.trigger_format);
input_file = sprintf('%s/%s/events-LaneChange.csv', output_dir, subject);
[data.sim.lane_change.values, data.sim.lane_change.hdr] = import_log(input_file, params.simlog.lanechange_format);
input_file = sprintf('%s/%s/events-IncreaseSpeed.csv', output_dir, subject);
[data.sim.accelerate.values, data.sim.accelerate.hdr] = import_log(input_file, params.simlog.accelerate_format);
input_file = sprintf('%s/%s/events-RewardAssessed.csv', output_dir, subject);
[data.sim.reward.values, data.sim.reward.hdr] = import_log(input_file, params.simlog.reward_format);
input_file = sprintf('%s/%s/events-RewardDisplayed.csv', output_dir, subject);
[data.sim.rewarddisp.values, data.sim.rewarddisp.hdr] = import_log(input_file, params.simlog.rewarddisp_format);
input_file = sprintf('%s/%s/events-SimulatorState.csv', output_dir, subject);
[data.sim.states.values, data.sim.states.hdr] = import_log(input_file, params.simlog.state_format);
input_file = sprintf('%s/%s/events-SimulatorEnded.csv', output_dir, subject);
[data.sim.simended.values, data.sim.simended.hdr] = import_log(input_file, params.simlog.simended_format);

% Synch time series
% Construct matrix with tracker time, linked to event distance and event
% simulation time.

log_ids = data.eye.log.messages{3}; %cell2mat(data.eye.log.messages(:,3));
log_times = data.eye.log.messages{1}; %cell2mat(data.eye.log.messages(:,1));
log_times = log_times / 1000;
log_times = log_times - results.t_start;

% Triggers (distance is screwed up on ver. 1.0 of ar.OpenDS)
idx_dist = find(strcmp(data.sim.triggers.hdr,'Distance'));
idx_type = find(strcmp(data.sim.triggers.hdr,'Action'));

% First instance of iterate cycle gives cycle total distance
types = data.sim.triggers.values{idx_type};
dists = data.sim.triggers.values{idx_dist};
if lane_dist < 0
    idx = find(strcmp(types,'iterate-cycle'));
    lane_dist = dists(idx(1));
end

% LaneChanges
idx_simtime = find(strcmp(data.sim.lane_change.hdr,'SimTime'));
idx_dist = find(strcmp(data.sim.lane_change.hdr,'Sim:DriverCar:LaneDistance'));
idx_cycle = find(strcmp(data.sim.lane_change.hdr,'Sim:Game:Cycle'));
idx_repeat = find(strcmp(data.sim.lane_change.hdr,'Sim:Game:Repeat'));
idx_id = find(strcmp(data.sim.lane_change.hdr,'LogId'));

ids = data.sim.lane_change.values{idx_id};
dists = data.sim.lane_change.values{idx_dist};
cycles = data.sim.lane_change.values{idx_cycle};
repeats = data.sim.lane_change.values{idx_repeat};
simtime = data.sim.lane_change.values{idx_simtime};

% IncreaseSpeed events
idx_simtime = find(strcmp(data.sim.accelerate.hdr,'SimTime'));
idx_dist = find(strcmp(data.sim.accelerate.hdr,'Sim:DriverCar:LaneDistance'));
idx_cycle = find(strcmp(data.sim.accelerate.hdr,'Sim:Game:Cycle'));
idx_repeat = find(strcmp(data.sim.accelerate.hdr,'Sim:Game:Repeat'));
idx_id = find(strcmp(data.sim.accelerate.hdr,'LogId'));

ids = [ids;data.sim.accelerate.values{idx_id}];
dists = [dists;data.sim.accelerate.values{idx_dist}];
cycles = [cycles;data.sim.accelerate.values{idx_cycle}];
repeats = [repeats;data.sim.accelerate.values{idx_repeat}];
simtime = [simtime;data.sim.accelerate.values{idx_simtime}];

% SimulatorState events
if ~isempty(data.sim.states.values{1})
    idx_simtime = find(strcmp(data.sim.states.hdr,'SimTime'));
    idx_dist = find(strcmp(data.sim.states.hdr,'Sim:DriverCar:LaneDistance'));
    idx_cycle = find(strcmp(data.sim.states.hdr,'Sim:Game:Cycle'));
    idx_repeat = find(strcmp(data.sim.states.hdr,'Sim:Game:Repeat'));
    idx_id = find(strcmp(data.sim.states.hdr,'LogId'));
    
    ids = [ids;data.sim.states.values{idx_id}];
    dists = [dists;data.sim.states.values{idx_dist}];
    cycles = [cycles;data.sim.states.values{idx_cycle}];
    repeats = [repeats;data.sim.states.values{idx_repeat}];
    simtime = [simtime;data.sim.states.values{idx_simtime}];
    
end

cycles(cycles==0)=1;
dists(dists<0)=0;

% Get AdjSerialByte
ab_ids = data.sim.events.values{ find(strcmp(data.sim.events.hdr, 'LogId')) };
adjbytes = data.sim.events.values{ find(strcmp(data.sim.events.hdr, 'AdjSerialByte')) };

% Match IDs to times
sim2track.hdr=[{'LogId'},{'TrackerTime'},{'Cycle'},{'Repeat'},{'LaneDist'},{'SimTime'},{'AdjSerialByte'}];
M = cell(length(ids),7);
pff = nan(length(ids),1);
isok = true(length(ids),1);
for j = 1 : length(ids)
    
    id = ids(j);
    idx = find(log_ids==id);
    
    idx2 = find(ab_ids == id);
    
    if isempty(idx) || isempty(id)
       isok(j) = false; 
    else
       try
           pff(j) = idx;
       catch err
           a = 0;
       end
       if isempty(idx2)
           M(j,:) = [{id},{log_times(idx)},{cycles(j)},{repeats(j)},{dists(j)},{simtime(j)},{0}];
       else
           M(j,:) = [{id},{log_times(idx)},{cycles(j)},{repeats(j)},{dists(j)},{simtime(j)},{adjbytes(idx2)}];
       end
    end
    
end

M = M(isok,:);
ids = ids(isok);

[~,idx] = sort(ids);
M = M(idx,:);

% Insert round 
j = 3; c = 1; r = 1;
last_c = 1; last_r = 1;
M2 = cell(0,7);
N = size(M,1);
while j <= N
    c = M{j,3};
    r = M{j,4};
    
    if c > last_c && r == 1
        k=j-1;
%         fprintf('k=%d | c=%d | r=%d\n', k, c, r);
        % This is a new cycle; interpolate to end
        if M{k,5} < lane_dist
           ddist1 = M{j,5} - M{k,5} + lane_dist;
           ddist2 = lane_dist - M{k,5};
           ratio = ddist2 / ddist1;
           dsimtime = M{j,6} - M{k,6};
           dsimtime = M{k,6} + dsimtime * ratio;
           dlogtime = M{j,2} - M{k,2};
           dlogtime = M{k,2} + dlogtime * ratio;
           
           % Insert
           M2(end+1,:) = [{M{k,1}+1},{dlogtime},{last_c},{last_r},{lane_dist},{dsimtime},{0}];
           M2(end+1,:) = [{M{k,1}+2},{dlogtime},{c},{r},{0},{dsimtime},{0}];
            
        end
        
    end
    
    last_c = c;
    last_r = r;
    j = j + 1;
end

% Append last part to end
j = N;
ddist1 = M{j,5} - M{j-1,5};
ddist2 = lane_dist - M{j,5};
dsimtime = M{j,6} - M{j-1,6};
dlogtime = M{j,2} - M{j-1,2};
dsimtime = ddist2 / ddist1 * dsimtime;
dlogtime = ddist2 / ddist1 * dlogtime;
M2(end+1,:) = [{M{j,1}+1},{M{j,2}+dlogtime},{last_c},{last_r},{lane_dist},{M{j,6}+dsimtime},{0}];

M = [M;M2];
[~,idx] = sort(cell2mat(M(:,1)));
M = M(idx,:);

idx = true(size(M,1),1);
% Remove double entries
for j = 2 : size(M,1)
    if abs(M{j,5}-M{j-1,5})<0.0001
        idx(j)=false;
    end
end

M = M(idx,:);

% Remove badly sorted times
while true
    idx = [];
    cycles = cell2mat(M(:,3));
    N_cycles = max(cycles);
    for cycle = 1:N_cycles
        idxc = find(cycles==cycle);
        S = M(idxc,:);
        repeats = cell2mat(S(:,4));
        N_repeats = max(repeats);
        for repeat = 1:N_repeats
            idxr = find(repeats==repeat);
            S2 = S(idxr,:);
            dists = cell2mat(S2(:,5));
            dd = diff(dists);
            idxd = find(dd<=0);
            if ~isempty(idxd)
               idx(end+1:end+length(idxd)) = idxc(1) + idxr(1) - 2 + idxd;

            end
        end
    end

    if isempty(idx), break; end

    N = size(M,1);
    idxx=true(N,1);
    idxx(idx)=false;
    M=M(idxx,:);
end

sim2track.matrix = M;

% Interpolate distance and sim-time to tracker time series
% Triggers (distance is screwed up on ver. 1.0 of ar.OpenDS)
idx_time = find(strcmp(data.sim.triggers.hdr,'SimTime'));
idx_type = find(strcmp(data.sim.triggers.hdr,'Action'));

% First instance of iterate cycle gives cycle total distance
types = data.sim.triggers.values{idx_type};
times = data.sim.triggers.values{idx_time};
idx = find(strcmp(types,'iterate-cycle'));
if n_rounds > 0
    if length(idx) == n_rounds - 1
        % Don't remove last index; simulation was likely terminated
        % prematurely
        warning('Using fixed number of rounds [%s] = %d', subject, n_rounds);
    elseif length(idx) == n_rounds
        idx = idx(1:end-1);
    else
        warning('Invalid number of rounds [%s] = %d, %d', subject, n_rounds, length(idx)-1)
    end
else
    idx = idx(1:end-1); % Last iteration is the end point
end

sim2track.cycle_times = interpolate_log_times( sim2track.matrix, sim2track.hdr, times(idx) );

% Rewards
idx_time = find(strcmp(data.sim.reward.hdr,'SimTime'));
idx_type = find(strcmp(data.sim.reward.hdr,'Type'));
idx_magnitude = find(strcmp(data.sim.reward.hdr,'Magnitude'));
types = data.sim.reward.values{idx_type};
times = data.sim.reward.values{idx_time};
magnitudes = data.sim.reward.values{idx_magnitude};
idx_collision = find(strcmp(types,'collision'));
sim2track.repeat_times = interpolate_log_times( sim2track.matrix, sim2track.hdr, times(idx_collision) );

coll_times = times(idx_collision);
coll_mags = magnitudes(idx_collision);
coll_types = types(idx_collision);

% Stupidness: find the nearest RewardDisplayed because something is wrong
% with rewards
idx_type = find(strcmp(data.sim.rewarddisp.hdr,'RewardType'));
idx_msg = find(strcmp(data.sim.rewarddisp.hdr,'Message'));
idx_id = find(strcmp(data.sim.rewarddisp.hdr,'RewardId'));
idx_time = find(strcmp(data.sim.rewarddisp.hdr,'SimTime'));
ids = data.sim.rewarddisp.values{idx_id};
magnitudes = data.sim.rewarddisp.values{idx_msg};
types = data.sim.rewarddisp.values{idx_type};
dtimes = data.sim.rewarddisp.values{idx_time};

mtimes = zeros(length(dtimes),1);

% map by times
for i = 1 : length(dtimes)
   idxt = find(times<=dtimes(i));
   mtimes(i) = times(idxt(end));
end

% Add collisions
mtimes = [mtimes;coll_times];
types = [types;coll_types];
magnitudes = [magnitudes;coll_mags];

idx = find(strcmp(types,'overtake'));
sim2track.overtake_times = interpolate_log_times( sim2track.matrix, sim2track.hdr, mtimes(idx) );

% All overtake-relevant rewards
% NB: Using sources here because "types" was incorrect in early version of
% OpenDS simulator
idx = [find(strcmp(types,'overtake')); ...
       find(strcmp(types,'overtake-traffic')); ...
       find(strcmp(types,'too-close')); ...
       find(strcmp(types,'collision'));];
sim2track.reward_times = interpolate_log_times( sim2track.matrix, sim2track.hdr, mtimes(idx) );

sim2track.reward_magnitudes = magnitudes(idx);
[sim2track.reward_times, idx] = sort(sim2track.reward_times);
sim2track.reward_magnitudes = sim2track.reward_magnitudes(idx);

% SimulationEnded event
idx_time = find(strcmp(data.sim.simended.hdr,'SimTime'));
sim2track.simended_time = interpolate_log_times( sim2track.matrix, sim2track.hdr, data.sim.simended.values{idx_time}, true );

% Lane changes
idx_dir = find(strcmp(data.sim.lane_change.hdr,'Direction'));
idx_simtime = find(strcmp(data.sim.lane_change.hdr,'SimTime'));
times = data.sim.lane_change.values{idx_simtime};
dirs = data.sim.lane_change.values{idx_dir};

idx = find(strcmp(dirs,'left'));
sim2track.left_change_times = interpolate_log_times( sim2track.matrix, sim2track.hdr, times(idx) );
sim2track.left_change_times = sim2track.left_change_times(~isnan(sim2track.left_change_times));
idx = find(strcmp(dirs,'right'));
sim2track.right_change_times = interpolate_log_times( sim2track.matrix, sim2track.hdr, times(idx) );
sim2track.right_change_times = sim2track.right_change_times(~isnan(sim2track.right_change_times));

% Baseline periods - assign to time series
sim2track.baseline = interpolate_baseline_intervals( params.baseline.intervals, ...
                                                     results.tgap, ...
                                                     sim2track.matrix, ...
                                                     sim2track.hdr );


results.sim2track = sim2track;


end

