function sim = get_simulation_events_eeg( outdir, preproc, proc, data )
%%%%%%%%%%%%%
% Loads simulation log events into a Matlab table and defines
% epochs.
%

% Load eye preprocessing results
preproc_results = load(sprintf('%s/results.mat', outdir), 'results');

% Load sequence difficulty ratings
seq_diff = load( proc.params.sequence_difficulty_file );

sim = [];

%Load lange change events
input_file = sprintf('%s/events-LaneChange.csv', outdir);
[values, hdr] = import_log(input_file, preproc.params.simlog.lanechange_format);

T = values{1};
mylog.times = double(T - data.eeg.t_start_sim); % Unit is ms
mylog.triggers = values{find(strcmp(hdr,'AdjSerialByte'))};
mylog.logid = values{find(strcmp(hdr,'LogId'))};
mylog.type = values{find(strcmp(hdr,'EventType'))};
mylog.lane_from = values{find(strcmp(hdr,'LaneFrom'))};
mylog.lane_to = values{find(strcmp(hdr,'LaneTo'))};
mylog.sim_time =values{find(strcmp(hdr,'SimulationTime'))};
mylog.cycle=values{find(strcmp(hdr,'Sim:Game:Cycle'))};
mylog.repeat=values{find(strcmp(hdr,'Sim:Game:Repeat'))};
mylog.lane_dist=values{find(strcmp(hdr,'Sim:DriverCar:LaneDistance'))};

mylog.difficulty = get_difficulty(mylog, seq_diff.difficulty);
mylog.outcomes = get_outcomes(mylog, preproc_results.results.sim2track);

sim.events = table(mylog.times, mylog.triggers, mylog.logid, mylog.type, mylog.lane_from, ...
                   mylog.lane_to, mylog.sim_time, mylog.cycle, mylog.repeat, mylog.lane_dist, ...
                   mylog.difficulty, mylog.outcomes);
sim.events.Properties.VariableNames = {'Time', 'Trigger', 'LogId', 'Type', 'LaneFrom', ...
                                       'LaneTo', 'SimTime', 'Cycle', 'Repeat', 'LaneDist', ...
                                       'Difficulty', 'Outcomes'};

% Define passing and baseline epochs
proc_results_file = sprintf('%s/processing_results.mat', outdir);
proc_results = load(proc_results_file);

sim.epochs.idx_baseline = proc_results.results.epochs.intervals.baseline.idx;
sim.epochs.idx_passing = proc_results.results.epochs.intervals.passing.idx;

min_epoch = 500;

sim.trl_epochs = zeros(0,5);
for b = 1 : size(sim.epochs.idx_baseline,1)
   idxb =  sim.epochs.idx_baseline(b,:);
   if idxb(2) - idxb(1) >= min_epoch
        sim.trl_epochs(end+1,:) = [idxb(1),idxb(2),0,1,1];
   end
end
for b = 1 : size(sim.epochs.idx_passing,1)
   idxb =  sim.epochs.idx_passing(b,:);
   if idxb(2) - idxb(1) >= min_epoch
        sim.trl_epochs(end+1,:) = [idxb(1),idxb(2),0,1,2];
   end
end



    function diffs = get_difficulty(mylog, difficulty)
        
        froms = cell2mat(difficulty(:,4));
        tos = cell2mat(difficulty(:,5));
        rounds = cell2mat(difficulty(:,1));
        diffs = nan(length(mylog.times),1);
       
        for i = 1 : length(mylog.times)
        
            dist = mylog.lane_dist(i) / proc.params.lane_dist;
            round = mylog.cycle(i);
            
            idx_round = find(rounds == round);
            
            if ~isempty(idx_round)
                idx_diff = find(dist > froms(idx_round) & dist < tos(idx_round),1);

                if ~isempty(idx_diff)
                    diffs(i) = difficulty{idx_round(idx_diff),6};
                end
            end
        end

    end

    function outcomes = get_outcomes(mylog, sim2track)
        
        maxdelta = 10000;
        outcomes = nan(length(mylog.times),1);
       
        for i = 1 : length(mylog.times)
        
            idx = find(sim2track.reward_times > mylog.times(i),1);
            if ~isempty(idx)
                i0 = max(1,idx-1); i1 = min(length(sim2track.reward_times),idx+1);
                idx2 = i0:i1;
                idxgt = abs(sim2track.reward_times(idx2) - mylog.times(i)) < maxdelta;
%                 if abs(sim2track.reward_times(idx) - mylog.times(i)) < maxdelta || ...
%                    (idx >1 && abs(sim2track.reward_times(idx-1) - mylog.times(i)) < maxdelta)
                if any(idxgt)
                    outcomes(i) = sim2track.reward_magnitudes(idx2(find(idxgt,1)));
                end
            else
                a=0;
            end
        end

    end


end