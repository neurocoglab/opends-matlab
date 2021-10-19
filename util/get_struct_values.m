function values = get_struct_values( my_struct, field, is_numeric )
    
    if nargin < 3
        is_numeric = false;
    end

    if is_numeric
        values = nan(length(my_struct),1);
        for s = 1 : length(my_struct)
            if ~isempty(my_struct(s).(field))
                values(s) = my_struct(s).(field);
            end
        end
    else
        values = cell(length(my_struct),1);
        for s = 1 : length(my_struct)
            values(s) = {my_struct(s).(field)};
        end
    end
end
