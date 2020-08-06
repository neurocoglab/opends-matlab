function [ result ] = interpolate_messagebuttons_sim ( data,button_map, M )
% interpolate_messagebuttons_sim 
%
% Interpolates message button presses from the matrix M into log times 
% and maps them to numeric responses

presses = data.sim.messagebutton.values.Source;
for i = 1 : length(presses)
    presses(i) = {strrep(presses{i},'ButtonPressDialog.','')};
end
times = data.sim.messagebutton.values.Millis - data.sim.t_start;

buttons = unique(button_map.Source);
result = [];
result.buttons = buttons;

for i = 1 : length(buttons)

    button = buttons{i};
    idx = find(strcmp(button_map.Source,button));
    button_map_i = button_map(idx,:);
    
    idx = find(strcmp(presses, button));
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
    result.(button).values = values(~isnan(values));
  
end

end

