function [ results ] = process_traffic_decision_events( data, results, params )
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
    idx1 = find(roadsigns.rounds==round_i);
    names_i = roadsigns.names(idx1);
    dirs = roadsigns.direction(idx1);
    idx2 = find(strcmp(names_i,name_i));
    dir_chosen(j) = {dir_i};
    dir_corr(j) = {lower(dirs{idx2})};
    correct(j) = strcmpi(dir_i, dir_corr(j));
    idx_c = find(events_c.rounds==round_i & events_c.repeats==repeat_i, 1);
    if ~isempty(idx_c)
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





end

