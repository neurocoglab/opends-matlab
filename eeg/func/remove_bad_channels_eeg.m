function [ data ] = remove_bad_channels_eeg( params, data, subject )

    % Find bad channels for this subject, if defined
    bad_channels = [];
    bad_channel_file = sprintf('%s/%s/%s', params.io.input_dir, ...
                                           params.io.metadata_dir, ...
                                           params.eeg.bad_channel_file);
    
    if exist(bad_channel_file, 'file')
       opts = detectImportOptions(bad_channel_file);
       opts = setvartype(opts, 'char');
       bad_channels = readtable(bad_channel_file, opts);
       bad_channels.Properties.VariableNames = [{'Subject'},{'Channel'}];
    end

    % Check that data
    data.eeg.bad_channels = [];

    % Remove bad channels
    if ~isempty(bad_channels)
        T = bad_channels{strcmp(bad_channels.Subject, subject),2};
        data.eeg.bad_channels = T;
        if ~isempty(T)
            cfg = [];
            cfg.channel = data.eeg.ft.label;
            idx_rem = [];
            for c = 1 : length(T)
                idx_rem = [idx_rem find(strcmp(cfg.channel,T{c}))];
            end
            cfg.channel(idx_rem) = [];
            [~,data.eeg.ft] = evalc('ft_selectdata(cfg, data.eeg.ft);');

%             fprintf('%s: Removed %d bad channels.\n', subject, length(T));
        end
    end
    
    % Interpolate over bad channels
    cfg = [];
    cfg.badchannel = data.eeg.bad_channels;
    cfg.layout = 'biosemi64.lay';
    cfg_nbr = [];
    cfg_nbr.method = 'distance';
    cfg_nbr.layout = 'biosemi64.lay';
    cfg_nbr.channel = data.eeg.eeg_channels;
    
    [~,cfg.neighbours] = evalc('ft_prepare_neighbours(cfg_nbr);');
    if ~isempty(cfg.badchannel)
        [~,data.eeg.ft] = evalc('ft_channelrepair(cfg, data.eeg.ft);');
    end

end

