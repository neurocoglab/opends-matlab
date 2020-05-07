function [ result ] = interpolate_messagebuttons_sim ( data, button_map, M )
% interpolate_messagebuttons_sim 
%
% Interpolates message button presses from the matrix M into log times 
% and maps them to numeric responses

idx_remove = [];

% Deal with repeated message log events [BUG]
T = data.sim.messagebutton.values;
for i = 2 : height(T)
    if strcmp(T.Message{i},T.Message{i-1}) && T.SimTime(i)==T.SimTime(i-1)
        idx_remove = [idx_remove i];
    end
end

T(idx_remove,:) = [];

presses = T.Source;
for i = 1 : length(presses)
    presses(i) = {strrep(presses{i},'ButtonPressDialog.','')};
    % For early iterations, strip the trailing integers
    [~,st] = str2num(presses{i}(end));
    if st
        str = presses{i};
        idx = strfind(str,'_');
        presses(i) = {str(1:idx(end)-1)};
    end
end

times = T.Millis - data.sim.t_start;

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
    result.(button).keys = T.KeyPressed(idx);
    % Omit empty responsess
    idx2 = idx2 & ~cellfun('isempty',result.(button).keys);
    idx = idx(idx2);
    result.(button).keys = result.(button).keys(idx2);
    result.(button).times = result.(button).times(idx2);
    result.(button).durations = T.Duration(idx);
    result.(button).messages = T.Message(idx);
    result.(button).rounds = T.Sim_Game_Cycle(idx);
    result.(button).repeats = T.Sim_Game_Repeat(idx);
    result.(button).serial_byte = T.AdjSerialByte(idx);
    values = nan(length(idx),1);
    for j = 1 : length( idx )
        key = T.KeyPressed(idx(j));
        if isempty(key{:})
            values(j) = -1;
        else
            try
            values(j) = button_map_i.Value(find(strcmp(button_map_i.Key,key),1));
            catch
               a=0; 
            end
        end
    end
    result.(button).values = values;
  
end

end

