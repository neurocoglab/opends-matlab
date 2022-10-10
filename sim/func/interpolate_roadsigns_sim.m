function [ result ] = interpolate_roadsigns_sim( data, roadsign_map, M, sim_version )
% interpolate_roadsigns_sim
%
% Interpolates road sign changes from the matrix M into log times 
% and maps them to left/right, first/second


result = [];
result.times = [];
result.direction = {};
result.order = [];
result.names = {};
result.rounds = [];
result.repeats = [];

sign_ids = data.sim.roadsigns.values.Source;

if height(sign_ids) == 0
    warning('No road sign data found');
   return; 
end


messages = data.sim.roadsigns.values.CurrentText;
rounds = data.sim.roadsigns.values.Sim_Game_Cycle;
repeats = data.sim.roadsigns.values.Sim_Game_Repeat;
times = data.sim.roadsigns.values.Millis - data.sim.t_start;

ids = roadsign_map.Id;


% Add all ids to list
for i = 1 : length(ids)
    id_i = ids{i};
    idx_i = strcmp(roadsign_map.Id,id_i);
    idx = find(strcmp(sign_ids, id_i));
    dir = roadsign_map.Direction(idx_i,:);
    order = roadsign_map.Position(idx_i,:);
    
    mm = messages(idx);
    for j = 1 : length(mm)
        msg = mm{j};
        idx_j = idx(j);
        parts = strsplit(msg, ' | ');
        for k = 1 : length(parts)
            result.times = [result.times;times(idx_j)];
            result.direction = [result.direction;dir];
            result.order = [result.order;order];
            result.names = [result.names;parts(k)];
            result.rounds = [result.rounds;rounds(idx_j)];
            result.repeats = [result.repeats;repeats(idx_j)];
        end
        
    end
    
%     result.times = [result.times;times(idx)];
%     result.direction = [result.direction;repmat(dir,length(idx),1)];
%     result.order = [result.order;repmat(order,length(idx),1)];
%     result.message = messages(idx);
    
end

[~,idx_s] = sort(result.times);
result.times = result.times(idx_s);
result.direction = result.direction(idx_s);
result.order = result.order(idx_s);
result.names = result.names(idx_s);
result.rounds = result.rounds(idx_s);
result.repeats = result.repeats(idx_s);

% times(times < 0) = 0; % Move events happening before t=0 to t=0

result.times = interpolate_log_times_sim( M, result.times, true );

end

