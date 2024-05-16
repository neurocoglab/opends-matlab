function [ T_cycles ] = get_cycle_indices_sim( params, data )
%GET_CYCLE_INDICES_SIM Get the start and stop trigger events (adjusted 
% serial bytes) of cycles and repeats for a simulation log
%
% Output is a table with trigger bytes for the start and stop of
%  each cycle/repeat

T_start = data.sim.simstart.values;
T_resumed = data.sim.resumed.values;
T_end = data.sim.simended.values;

hdr = {'Cycle','Repeat','SerialByte','Millis'};
T_cycles = cell2table({1,1,1,T_start.Millis(1)}, 'VariableNames', hdr);

c = 1;
for i = 1 : height(T_resumed)
    c_i = T_resumed.Game_Cycle(i);
    if c_i == c + 1
        T_cycles = [T_cycles; cell2table({c_i, ...
                                           T_resumed.Game_Repeat(i), ...
                                           T_resumed.AdjSerialByte(i), ...
                                           T_resumed.Millis(i)}, ...
                                           'VariableNames', hdr)];
        c = c_i;
    end
    
end

% End of last round is abs_max
if height(T_end) == 0
    warning('No SimulationEnded event found in log.');
else
    T_cycles = [T_cycles; cell2table({ T_cycles.Cycle(end), ...
                                       T_cycles.Repeat(end), ...
                                       T_end.AdjSerialByte(1), ...
                                       T_end.Millis(1)}, ...
                                       'VariableNames', hdr)];
end

end

