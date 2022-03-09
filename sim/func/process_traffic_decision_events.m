function [ results ] = process_traffic_decision_events( signals, data, results, params )
% Processes traffic decision button press events and associated confidence ratings
%

% Get times for traffic decision events and associated confidence events
events_d = data.sim.sim2track.messagebutton.direction_dialog;
events_c = data.sim.sim2track.messagebutton.confidence_dialog;
roadsigns = data.sim.sim2track.roadsigns;

N_ev = length(events_d.times);

% Deal with repeats in log...
cc = cell(N_ev,1);
for i = 1 : N_ev
    cc(i) = {sprintf('%d.%d',events_d.rounds(i),events_d.repeats(i))};
end
ucc = unique(cc);
idx=[];
for i = 1 : length(ucc)
    idx = [idx find(strcmp(cc,ucc{i}),1)];
end
idx = sort(idx);
N_ev = length(idx);

ids = repmat({data.subject},1,N_ev);

correct = false(N_ev,1);
dir_corr = cell(N_ev,1);
dir_chosen = cell(N_ev,1);
names = cell(N_ev,1);
confidence = ones(N_ev,1)*-1;
rounds = ones(N_ev,1)*-1;

% Compute right/wrong decision
for j = 1 : N_ev
    i = idx(j);
    round_i = events_d.rounds(i);
    rounds(j) = round_i;
    repeat_i = events_d.repeats(i);
    name_i = events_d.messages{i};
    name_i = name_i(14:end-1);
    names(j) = {name_i};
    dir_i = events_d.keys{i};
    dir_i = dir_i(8:end);
    
    % Get road signs matching this round
    idx1 = find(roadsigns.rounds==round_i);
    names_i = roadsigns.names(idx1);
    
    % Get correct decision direction (left/right)
    dirs = roadsigns.direction(idx1);
    idx2 = find(strcmp(names_i,name_i), 1);
    dir_chosen(j) = {dir_i};
    dir_corr(j) = {lower(dirs{idx2})};
    correct(j) = strcmpi(dir_i, dir_corr(j));
    idx_c = find(events_c.rounds==round_i & events_c.repeats==repeat_i, 1);
    if ~isempty(idx_c)
        % For debug
        if idx_c > length(events_c.values)
           a = 0; 
        end
        confidence(j) = events_c.values(idx_c);
    end
end

C = [ids(:), num2cell(events_d.times(idx)),names(:),dir_chosen(:),dir_corr(:), ...
     num2cell(correct(:)),num2cell(confidence(:)),num2cell(events_d.durations(idx)), ...
     num2cell(rounds(:))];

T = cell2table(C,'VariableNames',[{'Subject'},{'Time'},{'Name'},{'Direction'},{'Dir_corr'}, ...
                                  {'Correct'},{'Confidence'},{'RT'},{'Round'}]);

% Output as CSV table
results.sim.events.traffic_decision.values = T;

if params.sim.events.traffic_decision.write_to_file
    subj_dir = sprintf('%s/%s/%s', params.io.output_dir, ...
                                   data.subject, ...
                                   params.sim.sub_dir);
    if ~exist(subj_dir, 'dir'), mkdir(subj_dir); end
    output_file = sprintf('%s/traffic_decisions.csv', subj_dir);
   
    writetable(T, output_file);
    
end

% Get time-locked data - Direction dialog
events = T.Time;

[tlocked,tstart] = get_tlocked(signals.pdz, signals.t_pd, events, params.eye.events.traffic_decision);
% if params.eye.events.tlock_params.apply
%     tlock_params = get_tlocked_parameters(tlocked, params.eye.events.traffic_decision);
%     results.eye.events.traffic_decision.tlock_params = tlock_params;
% end
clear baseline_stats;
baseline_stats.mean = nan(length(tstart),1);
baseline_stats.std = nan(length(tstart),1);
tlocked_bl = zeros(size(tlocked));
for i = 1 : length(tstart)
   idx_bl = find(signals.t_baselines < tstart(i),1,'last');
   if ~isempty(idx_bl)
      xx = signals.pdz(signals.idx_baseline(idx_bl,1):signals.idx_baseline(idx_bl,2));
      baseline_stats.mean(i) = mean(xx);
      baseline_stats.std(i) = std(xx);
      tlocked_bl(i,:) = tlocked(i,:) - baseline_stats.mean(i);
   end
end
results.eye.events.traffic_decision.baseline_stats = baseline_stats;
results.eye.events.traffic_decision.events = events;
results.eye.events.traffic_decision.tlocked = tlocked;
results.eye.events.traffic_decision.tlocked_bl = tlocked_bl;
results.eye.events.traffic_decision.t = -params.eye.events.traffic_decision.prepost(1):params.eye.events.traffic_decision.prepost(2);
results.eye.events.traffic_decision.t = results.eye.events.traffic_decision.t/params.eye.Fs;
results.eye.events.traffic_decision.diffs = zeros(length(events),1);
results.eye.events.traffic_decision.outcomes = nan(length(events),1);

if params.eye.events.traffic_decision.confidence.apply
    results.eye.events.traffic_decision.confidence = nan(length(events),1);
end
if params.eye.events.traffic_decision.correct.apply
    results.eye.events.traffic_decision.correct = nan(length(events),1);
end

for i = 1 : length(events)
    if params.eye.events.traffic_decision.confidence.apply
        results.eye.events.traffic_decision.confidence(i) = (confidence(i)==3) + 1;
    end
    if params.eye.events.traffic_decision.correct.apply
        results.eye.events.traffic_decision.correct(i) = (correct(i)) + 1;
    end
end

% Get corresponding random time-locked events for statistical comparison
% N_event = length(events);
% N_rand = params.eye.events.random.N;
% events = get_random_events( events, N_rand, ...
%                             [results.eye.events.traffic_decision.t(1) ...
%                             results.eye.events.traffic_decision.t(end)], ...
%                             signals.idx_baseline );
% T = zeros(N_rand,N_event,length(results.eye.events.traffic_decision.t));
% T_params = [];
% T_params.slope = zeros(N_rand,N_event);
% T_params.amplitude = zeros(N_rand,N_event);
% for i = 1 : N_rand
%     tlocked = get_tlocked(signals.pdz, signals.t_pd, events(:,i), params.eye.events.traffic_decision);
% %     if params.eye.events.tlock_params.apply
% %         tlock_params = get_tlocked_parameters(tlocked, params.eye.events.traffic_decision);
% %         T_params.slope(i,:) = tlock_params.slopes;
% %         T_params.amplitude(i,:) = tlock_params.amplitudes;
% %     end
% end
% results.eye.events.traffic_decision.tlocked_bl2 = squeeze(nanmean(T,1));
% if params.eye.events.tlock_params.apply
%     results.eye.events.traffic_decision.tlocked_params_bl2 = T_params;
% end

end

