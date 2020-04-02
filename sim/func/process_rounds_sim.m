function [ data ] = process_rounds_sim( data, params )
% PROCESS_ROUNDS - Process the simulation rounds using sim events,
%                  mapping simulation time to eye tracker time. Also
%                  maps baseline distances to eye tracker time series.

% Load blink-processed data
% infile = sprintf('%s/saccades.mat',sim_dir);
% load(infile);

% Does this subject have a lane distance specified? (Fix for bug in early
% version)
lane_dist = -1;
n_rounds = -1;

if isfield(params.sim.rounds.fixdist, 'values') && height(params.sim.rounds.fixdist.values) > 0
    subjs = params.sim.rounds.fixdist.values.subject;
    idx = find(strcmp(subjs, data.subject));
    if ~isempty(idx)
       T = params.sim.rounds.fixdist.values(idx(1),:);
       lane_dist = T.distance;
       n_rounds = T.rounds;
    end
end
 
sim_dir = sprintf('%s/%s/%s', params.io.output_dir, data.subject, params.sim.sub_dir);

% Load TriggerActuated | LaneChange | IncreaseSpeed events from simulation log
input_file = sprintf('%s/events-All.csv', sim_dir);
data.sim.events.values = import_log_sim(input_file, params.sim.log.events_format);
input_file = sprintf('%s/events-SimulatorStarted.csv', sim_dir);
data.sim.simstart.values = import_log_sim(input_file, params.sim.log.start_format);
input_file = sprintf('%s/events-TriggerActuated.csv', sim_dir);
data.sim.triggers.values = import_log_sim(input_file, params.sim.log.trigger_format);
input_file = sprintf('%s/events-LaneChange.csv', sim_dir);
data.sim.lane_change.values = import_log_sim(input_file, params.sim.log.lanechange_format);
input_file = sprintf('%s/events-IncreaseSpeed.csv', sim_dir);
data.sim.accelerate.values = import_log_sim(input_file, params.sim.log.accelerate_format);
input_file = sprintf('%s/events-RewardAssessed.csv', sim_dir);
data.sim.reward.values = import_log_sim(input_file, params.sim.log.reward_format);
input_file = sprintf('%s/events-RewardDisplayed.csv', sim_dir);
data.sim.rewarddisp.values = import_log_sim(input_file, params.sim.log.rewarddisp_format);
input_file = sprintf('%s/events-SimulatorState.csv', sim_dir);
data.sim.states.values = import_log_sim(input_file, params.sim.log.state_format);
input_file = sprintf('%s/events-SimulatorEnded.csv', sim_dir);
data.sim.simended.values = import_log_sim(input_file, params.sim.log.simended_format);
if params.sim.rounds.messagebutton.apply
    input_file = sprintf('%s/events-MessageButtonPress.csv', sim_dir);
    data.sim.messagebutton.values = import_log_sim(input_file, params.sim.log.messagebutton_format);
end

if params.sim.rounds.roadsign.apply
    input_file = sprintf('%s/events-SignTextChanged.csv', sim_dir);
    data.sim.roadsigns.values = import_log_sim(input_file, params.sim.log.roadsign_format);
end

% Synch time series
% Construct matrix with tracker time, linked to event distance and event
% simulation time.

log_ids = data.eye.log.messages.LogId; 
log_times = data.eye.log.messages.Time;

% Zero log times on first trigger
log_times = log_times - data.eye.t_start_sim;

% Triggers (distance is screwed up on ver. 1.0 of ar.OpenDS)

% First instance of iterate cycle gives cycle total distance
types = data.sim.triggers.values.Action;
dists = data.sim.triggers.values.Distance;

if lane_dist < 0
    idx = [find(strcmp(types,'iterate-cycle'));find(strcmp(types,'iterate-cycle-reset'))];
    sort(idx);
    lane_dist = dists(idx(1));
end

% Define simulation time offset (start to first trigger) in seconds
% Note: simulation started trigger typically is not inserted in time
% series, so we have to synchronise to the first button trigger event

% LaneChanges
ids = data.sim.lane_change.values.LogId; 
dists = data.sim.lane_change.values.Sim_DriverCar_LaneDistance; 
cycles = data.sim.lane_change.values.Sim_Game_Cycle;
repeats = data.sim.lane_change.values.Sim_Game_Repeat; 
simtime = data.sim.lane_change.values.Millis;

% IncreaseSpeed events
ids = [ids;data.sim.accelerate.values.LogId]; 
dists = [dists;data.sim.accelerate.values.Sim_DriverCar_LaneDistance]; 
cycles = [cycles;data.sim.accelerate.values.Sim_Game_Cycle];
repeats = [repeats;data.sim.accelerate.values.Sim_Game_Repeat]; 
simtime = [simtime;data.sim.accelerate.values.Millis];

% SimulatorState events
if height(data.sim.states.values) > 0
    ids = [ids;data.sim.states.values.LogId]; 
    dists = [dists;data.sim.states.values.Sim_DriverCar_LaneDistance]; 
    cycles = [cycles;data.sim.states.values.Sim_Game_Cycle];
    repeats = [repeats;data.sim.states.values.Sim_Game_Repeat]; 
    simtime = [simtime;data.sim.states.values.Millis];
end

simtime = simtime - data.sim.t_start;

cycles(cycles==0)=1;
dists(dists<0)=0;

% Get AdjSerialByte
ab_ids = data.sim.events.values.LogId;
adjbytes = data.sim.events.values.AdjSerialByte;

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
j = 3;
last_c = 1; last_r = 1;
M2 = cell(0,7);
N = size(M,1);
while j <= N
    c = M{j,3};
    r = M{j,4};
    
    if c > last_c && r == 1
        k=j-1;
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

sim2track.matrix = cell2table(M, 'VariableNames', sim2track.hdr);

if params.sim.rounds.fix_repeats
   sim2track.matrix = fix_repeats( sim2track.matrix );
end

% Interpolate distance and sim-time to tracker time series
% Triggers (distance is screwed up on ver. 1.0 of ar.OpenDS)

% First instance of iterate cycle gives cycle total distance
types = data.sim.triggers.values.Action;
times = data.sim.triggers.values.Millis - data.sim.t_start;

idx = find(strcmp(types,'iterate-cycle') | strcmp(types,'iterate-cycle-reset'));
if n_rounds > 0
    if length(idx) == n_rounds - 1
        % Don't remove last index; simulation was likely terminated
        % prematurely
        warning('Using fixed number of rounds [%s] = %d', data.subject, n_rounds);
    elseif length(idx) == n_rounds
        idx = idx(1:end-1);
    else
        warning('Invalid number of rounds [%s] = %d, %d', data.subject, n_rounds, length(idx)-1)
    end
else
    idx = idx(1:end-1); % Last iteration is the end point
end

sim2track.cycle_times = interpolate_log_times_sim( sim2track.matrix, times(idx) );

% Rewards
types = data.sim.reward.values.Type;
times = data.sim.reward.values.Millis - data.sim.t_start;
magnitudes = data.sim.reward.values.Magnitude;

idx_collision = find(strcmp(types,'collision'));
sim2track.repeat_times = interpolate_log_times_sim( sim2track.matrix, times(idx_collision) );

coll_times = times(idx_collision);
coll_mags = magnitudes(idx_collision);
coll_types = types(idx_collision);

% Stupidness: find the nearest RewardDisplayed because something is wrong
% with rewards
ids = data.sim.rewarddisp.values.RewardId;
magnitudes = data.sim.rewarddisp.values.Message;
types = data.sim.rewarddisp.values.RewardType;
dtimes = data.sim.rewarddisp.values.Millis - data.sim.t_start;

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
sim2track.overtake_times = interpolate_log_times_sim( sim2track.matrix, mtimes(idx) );

% All overtake-relevant rewards
% NB: Using sources here because "types" was incorrect in early version of
% OpenDS simulator
idx = [find(strcmp(types,'overtake')); ...
       find(strcmp(types,'overtake-traffic')); ...
       find(strcmp(types,'too-close')); ...
       find(strcmp(types,'collision'));];
sim2track.reward_times = interpolate_log_times_sim( sim2track.matrix, mtimes(idx) );

sim2track.reward_magnitudes = magnitudes(idx);
[sim2track.reward_times, idx] = sort(sim2track.reward_times);
sim2track.reward_magnitudes = sim2track.reward_magnitudes(idx);

% SimulationEnded event
sim2track.simended_time = interpolate_log_times_sim( sim2track.matrix, data.sim.simended.values.Millis - data.sim.t_start, true );

% Lane changes
times = data.sim.lane_change.values.Millis - data.sim.t_start;
dirs = data.sim.lane_change.values.Direction;

idx = find(strcmp(dirs,'left'));
sim2track.left_change_times = interpolate_log_times_sim( sim2track.matrix, times(idx) );
sim2track.left_change_times = sim2track.left_change_times(~isnan(sim2track.left_change_times));
idx = find(strcmp(dirs,'right'));
sim2track.right_change_times = interpolate_log_times_sim( sim2track.matrix, times(idx) );
sim2track.right_change_times = sim2track.right_change_times(~isnan(sim2track.right_change_times));

% Message button presses
if params.sim.rounds.messagebutton.apply
    message_button_map_file = sprintf('%s/%s/%s', params.io.input_dir, params.io.metadata_dir, params.sim.rounds.messagebutton.maps_file );
    message_button_map = readtable( message_button_map_file );
    sim2track.messagebutton = interpolate_messagebuttons_sim( data, ...
                                                              message_button_map, ...
                                                              sim2track.matrix ); 
end

% Road sign changes
if params.sim.rounds.roadsign.apply
    roadsign_map_file = sprintf('%s/%s/%s', params.io.input_dir, params.io.metadata_dir, params.sim.rounds.roadsign.maps_file );
    roadsign_map = readtable( roadsign_map_file );
    sim2track.roadsigns = interpolate_roadsigns_sim( data, ...
                                                     roadsign_map, ...
                                                     sim2track.matrix ); 
end

% Baseline periods - assign to time series
baseline_intervals_file = sprintf('%s/%s/%s', params.io.input_dir, params.io.metadata_dir, params.sim.baseline.intervals_file);
baseline_intervals = readtable(baseline_intervals_file);
sim2track.baseline = interpolate_baseline_intervals_sim( baseline_intervals, ...
                                                         data.eye.tgap, ...
                                                         sim2track.matrix );


data.sim.sim2track = sim2track;


    function Mnew = fix_repeats( matrix )
       
        hdr = matrix.Properties.VariableNames;
        Mnew = cell2table(cell(0,length(hdr)),'VariableNames',hdr);
        
        % First get rid of stupid rows
        for ii = 1 : height(matrix)-1
            
            last_dist = matrix.LaneDist(ii);
            last_cycle = matrix.Cycle(ii);
            this_dist = matrix.LaneDist(ii+1);
            this_cycle = matrix.Cycle(ii+1);
            
            if last_dist-this_dist > 100 && last_cycle == this_cycle
                aa=0;
            else
                Mnew = [Mnew;matrix(ii,:)];
            end
            
        end
        
        matrix = Mnew;
        rep = 1;
        
        for ii = 2 : height(matrix)
            
            last_dist = matrix.LaneDist(ii-1);
            last_cycle = matrix.Cycle(ii-1);
            this_dist = matrix.LaneDist(ii);
            this_cycle = matrix.Cycle(ii);
            
            if this_dist < last_dist
               if this_cycle == last_cycle
                   % this is a repeat
                   rep = rep + 1;
               else
                   % this is a new cycle
                   rep = 1;
               end
            end
            
            Mnew.Repeat(ii)=int32(rep);
            
        end
        
        
    end


end

