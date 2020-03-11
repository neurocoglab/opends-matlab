function [ data ] = remove_artifacts_eeg( params, data )
%REMOVE_ARTIFACTS_EEG Remove saccade and z-value artifacts from the EEG
% time series in "data"
%   

cfg = [];

% Artifact rejection based on ET-defined saccades and eyeblinks
if params.eeg.artifacts.eye.apply
   
    % Align time series in seconds
   t_eye = data.eye.t / 1000;
   t_eeg = data.eeg.ft.time{1};

   artifacts = zeros(0,2);

   % Blinks
   for i = 1 : length(data.eye.blinks.blink_ints)
        ints_i = data.eye.blinks.blink_ints{i};
        for j = 1 : size(ints_i,1)
            xs = [t_eye(ints_i(j,1)) t_eye(ints_i(j,1)+ints_i(j,2))];
            idx0 = findInSorted(t_eeg, xs(1));
            idx1 = findInSorted(t_eeg, xs(2));
            artifacts = [artifacts;[idx0 idx1]];
        end
   end 

   %Saccades
   idx_sacc = data.eye.saccades.saccades(:,1:2);
   for i = 1 : size(idx_sacc,1)
        xs = t_eye(idx_sacc(i,:));
        idx0 = findInSorted(t_eeg, xs(1));
        idx1 = findInSorted(t_eeg, xs(2));
        artifacts = [artifacts;[idx0 idx1]];
   end
   
   cfg.artfctdef.eog.artifact = artifacts;
   
end

if params.eeg.artifacts.zscore.apply
   
   cfgz = params.eeg.artifacts.zscore.cfg;
   cfgz.artfctdef.zvalue.channel = data.eeg.eeg_channels;
   
   [~, ~, artifacts] = evalc('ft_artifact_zvalue(cfgz, data.eeg.ft);');
   
   cfg.artfctdef.zvalue.artifact = artifacts;
    
end

if isempty(cfg)
   warning('\nNo artifact rejection specified.');
   return;
end

% Remove identified artifacts
cfg.artfctdef.reject          = 'nan';
cfg.artfctdef.minaccepttim    = params.eeg.artifacts.minaccepttim;
data.eeg.artifact.cfg         = cfg;
[~, data.eeg.ft] = evalc('ft_rejectartifact(cfg, data.eeg.ft);');

% Interpolate over artifacts
if params.eeg.artifacts.interpolate.apply
    cfg = params.eeg.artifacts.interpolate.cfg;
    cfg.feedback = 'no';
    [~, data.eeg.ft] = evalc('ft_interpolatenan_corrected(cfg, data.eeg.ft);');
end


end


