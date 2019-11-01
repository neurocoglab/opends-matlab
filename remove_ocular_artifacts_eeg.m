function [data, et_results] = remove_ocular_artifacts_eeg ( data, preproc, subject, show_plots )
%%%%%%%%%%%%%%%%%%%
% Loads previously processed eye data (see 'prerocessing.m' and 'processing.m') 
% for the given subject and removes ocular artifacts from the EEG data in 'data'.
%
%
%

if nargin < 4
   show_plots = false; 
end

outdir = sprintf('%s/%s/%s', preproc.params.root_dir, preproc.params.output_dir, subject);
figdir = sprintf('%s/figures', outdir);

%% Load eye data 
et_results_file = sprintf('%s/%s/%s/results.mat',preproc.params.root_dir, ...
                          preproc.params.output_dir, subject);
et_results = load(et_results_file); 

dt_eye = single(et_results.data.eye.log.messages{1}(1)/1000 - et_results.data.eye.t_start);
% Align data to first log events
et_results.results.t = et_results.results.t - dt_eye;

% cfg = data.eeg.cfg;
% cfg.continuous = 'yes';
% cfg.hpfilter = 'yes';
% cfg.hpfreq = 0.5;
% cfg.hpfiltord = 6;
% cfg.lpfilter = 'yes';
% cfg.lpfreq = 60;
% cfg.lpfiltord = 6;
% cfg.demean = 'no';


% Artifact rejection based on ET-defined saccades and eyeblinks
t_eye = et_results.results.t / 1000;
t_eeg = data.eeg.ft.time{1};

artifacts = zeros(0,2);

% Blinks
for i = 1 : length(et_results.results.blink_ints)
    ints_i = et_results.results.blink_ints{i};
    for j = 1 : size(ints_i,1)
        xs = [t_eye(ints_i(j,1)) t_eye(ints_i(j,1)+ints_i(j,2))];
        idx0 = findInSorted(t_eeg, xs(1));
        idx1 = findInSorted(t_eeg, xs(2));
        artifacts = [artifacts;[idx0 idx1]];
    end
end 

%Saccades
idx_sacc = et_results.results.saccades.saccades(:,1:2);
for i = 1158%1 : size(idx_sacc,1)
    xs = t_eye(idx_sacc(i,:));
    idx0 = findInSorted(t_eeg, xs(1));
    idx1 = findInSorted(t_eeg, xs(2));
    artifacts = [artifacts;[idx0 idx1]];
end

cfg_art = [];
cfg_art.continuous = 'yes';
cfg_art.artfctdef.zvalue.channel = data.eeg.eeg_channels;
cfg_art.artfctdef.zvalue.cutoff = 20;
cfg_art.artfctdef.zvalue.trlpadding = 0;
cfg_art.artfctdef.zvalue.fltpadding = 0;
cfg_art.artfctdef.zvalue.artpadding = 0.05;

% algorithmic parameters
%cfg_art.artfctdef.zvalue.cumulative = 'yes';
%cfg_art.artfctdef.zvalue.medianfilter = 'yes';
%cfg_art.artfctdef.zvalue.medianfiltord = 9;
cfg_art.artfctdef.zvalue.absdiff = 'yes'; % this is important!!!
%cfg_art.artfctdef.zvalue.detrend = 'yes';
cfg_art.artfctdef.zvalue.interactive = 'no';

%[~, artifactz] = ft_artifact_zvalue(cfg_art, data.eeg.ft);


[~, ~, artifactz] = evalc('ft_artifact_zvalue(cfg_art, data.eeg.ft);');
% Only accept 

cfg_art = [];
cfg_art.artfctdef.reject          = 'nan'; %params.eeg.artifacts.reject;
cfg_art.artfctdef.eog.artifact    = artifacts;
cfg_art.artfctdef.zvalue.artifact = artifactz;
%    cfg_art.artfctdef.minaccepttim    = params.eeg.artifacts.minaccepttim_fft;
data.eeg.cfg.artfctdef = cfg_art.artfctdef;
[~,data.eeg.ft] = evalc('ft_rejectartifact(cfg_art, data.eeg.ft);');

data.eeg.artifacts = artifacts;


if show_plots
    plot_eeg_blinks(et_results.results, data.eeg.ft, [{'vEOGover'},{'Fz'},{'FC5'},{'POz'},{'T7'}], 6, ...
       sprintf('%s/eeg_on_blinks.png', figdir) );
end



end