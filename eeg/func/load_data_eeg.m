function [ data ] = load_data_eeg( params, subject )

data = [];

if strcmp(params.eeg.convert.format,'brainvision')
    data = load_data_brainvision_eeg( params, subject );
elseif strcmp(params.eeg.convert.format,'biosemi')
    data = load_data_biosemi_eeg( params, subject );
else
    warning('Format %s not supported!', params.eeg.format);
end

end

