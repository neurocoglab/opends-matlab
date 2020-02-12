function [ data ] = remove_bad_channels_eeg( params, data )

    % Find bad channels for this subject, if defined
    bad_channels = [];

    if exist(params.eeg.bad_channel_file, 'file')
       opts = detectImportOptions(params.eeg.bad_channel_file);
       opts = setvartype(opts, 'char');
       bad_channels = readtable(params.eeg.bad_channel_file, opts);
       bad_channels.Properties.VariableNames = [{'Subject'},{'Channel'}];
    end

    % Check that data
    data.eeg.badchannels = [];

    % Remove bad channels
    if ~isempty(bad_channels)
        T = bad_channels{strcmp(bad_channels.Subject, subject),2};
        data.eeg.badchannels = T;
        if ~isempty(T)
            cfg = [];
            cfg.channel = data.eeg.ft.label;
            idx_rem = [];
            for c = 1 : length(T)
                idx_rem = [idx_rem find(strcmp(cfg.channel,T{c}))];
            end
            cfg.channel(idx_rem) = [];
            [~,data.eeg.ft] = evalc('ft_selectdata(cfg, data.eeg.ft)');

            fprintf('%s: Removed %d bad channels.\n', subject, length(T));
        end
    end
    
    % Interpolate over bad channels
    

end

