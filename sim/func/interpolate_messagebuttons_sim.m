function [ result ] = interpolate_messagebuttons_sim ( data, button_map, M, sim_version )
% interpolate_messagebuttons_sim 
%
% Interpolates message button presses from the matrix M into log times 
% and maps them to numeric responses

presses = data.sim.messagebutton.values.Source;
for i = 1 : length(presses)
    press_i = strrep(presses{i},'ButtonPressDialog.','');
    if compare_versions(sim_version, '1.8.3') < 0
        % Strip the trailing index for this version
        idx = strfind(press_i, '_');
        if ~isempty(idx)
            press_i = press_i(1:idx(end)-1);
        end
    end
    presses(i) = {press_i};
end
times = data.sim.messagebutton.values.Millis - data.sim.t_start;


[buttons,idx_b] = unique(button_map.Source);
if any(strcmp(button_map.Properties.VariableNames, 'Name'))
    button_names = button_map.Name(idx_b);
else
    button_names = buttons;
end
result = [];
result.buttons = buttons;

for i = 1 : length(button_names)

    button = buttons{i};
    button_name = button_names{i};
    idx = find(strcmp(button_map.Source,button));
    button_map_i = button_map(idx,:);
    
    idx = find(strcmp(presses, button_name));
    result.(button).times = interpolate_log_times_sim( M, times(idx) );
    idx2 = ~isnan(result.(button).times);
    idx = idx(idx2);
    result.(button).times = result.(button).times(idx2);
    result.(button).durations = data.sim.messagebutton.values.Duration(idx);
    result.(button).keys = data.sim.messagebutton.values.KeyPressed(idx);
    result.(button).keys = data.sim.messagebutton.values.KeyPressed(idx);
    result.(button).messages = data.sim.messagebutton.values.Message(idx);
    result.(button).rounds = data.sim.messagebutton.values.Sim_Game_Cycle(idx);
    result.(button).repeats = data.sim.messagebutton.values.Sim_Game_Repeat(idx);
    result.(button).serial_byte = data.sim.messagebutton.values.AdjSerialByte(idx);
    values = nan(length(idx),1);
    for j = 1 : length( idx )
        key = data.sim.messagebutton.values.KeyPressed(idx(j));
        idx2 = find(strcmp(button_map_i.Key,key),1);
        if isempty(idx2)
           continue; 
        end
        values(j) = button_map_i.Value(idx2);
    end
    
    % Remove rows where value is NaN
    idx = ~isnan(values);
    result.(button).times = result.(button).times(idx);
    result.(button).durations =  result.(button).durations(idx);
    result.(button).keys = result.(button).keys(idx);
    result.(button).messages = result.(button).messages(idx);
    result.(button).rounds = result.(button).rounds(idx);
    result.(button).repeats = result.(button).repeats(idx);
    result.(button).serial_byte = result.(button).serial_byte(idx);
    result.(button).values = values(idx);
  
end

end

