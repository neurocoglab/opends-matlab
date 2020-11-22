function [ T_cycles ] = get_cycle_triggers_sim( data )
%GET_CYCLE_INDICES_SIM Get the start and stop trigger events (adjusted 
% serial bytes) of cycles and repeats for a simulation log
%
% Output is a table with trigger bytes for the start and stop of
%  each cycle/repeat

T_triggers = data.sim.triggers.values;
T_resumed = data.sim.resumed.values;
T_collision = data.sim.reward.values;
T_collision = T_collision(strcmp(T_collision.Type, 'collision'),:);
abs_max = max(data.sim.events.values.AdjSerialByte);

% First round is defined by (SimulatorStart -> iterate-cycle OR 
%                            iterate-cycle-reset OR
%                            collision)
% Subsequent rounds defined by (SimulatorResumed [1st after iterate-cycle] ->
% iterate-cycle)

hdr = {'Cycle','Repeat','Position','AdjSerialByte'};
T_cycles = array2table({1,1,'start',1});
T_cycles.Properties.VariableNames = hdr;

idx_iter = find(strcmp(T_triggers.Action, 'iterate-cycle') | ...
                strcmp(T_triggers.Action, 'iterate-cycle-reset'));

asbs = T_triggers.AdjSerialByte(idx_iter);
iscollision = false(length(asbs),1);
if height(T_collision) > 0
    asbs = [asbs; T_collision.AdjSerialByte];
    iscollision = [iscollision; true(height(T_collision),1)];
end

[asbs, idx] = sort(asbs);
iscollision = iscollision(idx);
            
c = 1; r = 1;
for i = 1 : length(asbs)-1
    asb1 = asbs(i);
    T_rows = array2table({c,r,'stop',asb1});
    
    idx_r = find(T_resumed.AdjSerialByte > asb1, 1);
    asb2 = T_resumed.AdjSerialByte(idx_r);
    
    if iscollision(i)
        % This is a collision
        r = r + 1;
    else
        % No collision
        c = T_resumed.Game_Cycle(idx_r);
        r = T_resumed.Game_Repeat(idx_r);
    end
    
    T_rows = [T_rows; array2table({c,r,'start',asb2})];
    T_rows.Properties.VariableNames = hdr;
    T_cycles = [T_cycles; T_rows];
    
end

% End of last round is abs_max
T_rows = array2table({c,r,'stop',abs_max});
T_rows.Properties.VariableNames = hdr;
T_cycles = [T_cycles; T_rows];

end

