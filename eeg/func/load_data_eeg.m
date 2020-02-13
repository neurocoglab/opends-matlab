function [ data ] = load_data_eeg( params, subject )

data = [];

if strcmp(params.eeg.format,'brainvision')
    data = load_data_brainvision_eeg( params, subject );
else
    warning('Format %s not supported!', params.eeg.format);
end

end

