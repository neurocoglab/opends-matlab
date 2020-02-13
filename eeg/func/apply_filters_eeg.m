function [ data ] = apply_filters_eeg( params, data )
%APPLY_FILTERS_EEG 
%
% Applies sequential bandpass and/or notch filters to the
% data. For details see the Fieldtrip command "ft_preprocessing".
%
%
%

if params.eeg.bandpass.apply
   cfg = params.eeg.bandpass.cfg;
   cfg.channel = data.eeg.eeg_channels;
   cfg.bpfilter = 'no';
   cfg.bsfilter = 'no';
   [~,data.eeg.ft] = evalc('ft_preprocessing(cfg, data.eeg.ft);'); 
end

if params.eeg.notch.apply
   cfg = params.eeg.notch.cfg;
   cfg.channel = data.eeg.eeg_channels;
   cfg.hpfilter = 'no';
   cfg.lpfilter = 'no';
   cfg.bpfilter = 'no';
   cfg.bsfilter = 'yes';
   [~,data.eeg.ft] = evalc('ft_preprocessing(cfg, data.eeg.ft);'); 
end


end

