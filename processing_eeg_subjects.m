%% Perform some sanity checks on all subjects

clear;
load processing_eeg_params.mat
preproc = load('preproc_params_hd.mat');
proc = load('processing_params.mat');

addpath('../lib/fieldtrip-20190419/');
addpath('../lib/cPCOH/functions');
addpath('../lib/textprogressbar');

bad_channels = [];

show_plots = false;

if exist(params.eeg.bad_channel_file, 'file')
   opts = detectImportOptions(params.eeg.bad_channel_file);
   opts = setvartype(opts, 'char');
   bad_channels = readtable(params.eeg.bad_channel_file, opts);
   bad_channels.Properties.VariableNames = [{'Subject'},{'Channel'}];
end

load(proc.params.qc.file);
qc_score = cell2mat(qc_eeg(:,1));
idx = qc_score>1; %=proc.params.qc.cutoff;
subjects = qc_eeg(idx,2);

fig_output_dir = sprintf('%s/subjects', proc.params.output_dir);
if ~exist(fig_output_dir, 'dir')
   mkdir(fig_output_dir); 
end

channel = 'Fz';

if false
    
%% Plot time series - after ICA

close all;

h = figure;
h.Color = 'w';
p = 1;


for i = 1 : length(subjects)
    
    subplot(5,ceil(length(subjects)/5), p);
    p = p + 1;
    
    subject = subjects{i};
    
    outdir = sprintf('%s/%s/%s', preproc.params.root_dir, preproc.params.output_dir, subject);
    figdir = sprintf('%s/figures', outdir);
    
    data = load_eeg_data( params, preproc, subject );
    results_file = sprintf('%s/processing_results_eeg_ica.mat',outdir);
    T = load(results_file, 'data_flt');
    data.eeg.ft = T.data_flt;
    clear T;
    
    idx_c = find(strcmp(data.eeg.ft.label,'Fz'));
    X = squeeze(data.eeg.ft.trial{1}(idx_c,:));
    t = data.eeg.ft.time{1};
    idx_nan = isnan(X);
%     X2=X;
%     X2(:) = nan;
%     X2(idx_nan) = 0;
%     h2 = plot(X2);
%     h2.Color = 'r';
%     hold on;
    X(idx_nan) = zscore(X(idx_nan));
    X(abs(X)>50)=nan;
    X(idx_nan) = zscore(X(idx_nan));
    hh = plot(t,X);
    hh.Color = 'b';
    ax = gca;
    ax.XTick = [];
    ax.XTickLabel = [];

    ylim([-50 50]);
    xlim([min(t) max(t)]);
    
    title(subject);
    
    fprintf('Added subject %s\n', subject);
    
    
end

suptitle(sprintf('All time series after ICA - channel %s', channel));

%% Interactive artifact rejection

p = 1;

for i = 2%1 : length(subjects)
    
    subject = subjects{i};
    
    data = load_eeg_data( params, preproc, subject );
    results_file = sprintf('%s/processing_results_eeg_ica.mat',outdir);
    T = load(results_file, 'data_flt');
    data.eeg.ft = T.data_flt;
    clear T;
    
    
    cfg=[];
    cfg.continuous = 'yes';
    cfg.artfctdef.zvalue.channel = data.eeg.eeg_channels;
    cfg.artfctdef.zvalue.cutoff = 20;
    cfg.artfctdef.zvalue.trlpadding = 0;
    cfg.artfctdef.zvalue.artpadding = 0.05;
    cfg.artfctdef.zvalue.fltpadding = 0;

    % algorithmic parameters
    cfg.artfctdef.zvalue.cumulative = 'yes';
    cfg.artfctdef.zvalue.medianfilter = 'yes';
    cfg.artfctdef.zvalue.medianfiltord = 9;
    cfg.artfctdef.zvalue.absdiff = 'yes';
    cfg.artfctdef.zvalue.detrend = 'yes';

    % make the process interactive
    cfg.artfctdef.zvalue.interactive = 'yes';

    [cfg, artifactz] = ft_artifact_zvalue(cfg, data.eeg.ft);
    cfg_art = [];
    cfg_art.artfctdef.reject          = 'nan'; %params.eeg.artifacts.reject;
    cfg_art.artfctdef.zvalue.artifact = artifactz;
    
    data.eeg.ft2 = ft_rejectartifact(cfg_art, data.eeg.ft);
    
    h=figure;
    h.Color='w';
    
    subplot(5,ceil(length(subjects)/5), p); hold on
    p = p + 1;
    figure; hold on 
    idx_c = find(strcmp(data.eeg.ft.label,'Fz'));
    
    t = data.eeg.ft.time{1};
    plot(data.eeg.ft.time{1},data.eeg.ft.trial{1}(idx_c,:))
    plot(data.eeg.ft.time{1},data.eeg.ft2.trial{1}(idx_c,:))
    ylim([-50 50]);
    xlim([min(t) max(t)]);
    
    title(subject);
    
end

suptitle(sprintf('All time series after Z art rej - channel %s', channel));

%% Plot time series - after artifact rejection

close all;

h = figure;
h.Color = 'w';
p = 1;
channel = 'Fz';

for i = 1 : length(subjects)
    
    subplot(5,ceil(length(subjects)/5), p);
    p = p + 1;
    
    subject = subjects{i};
    
    outdir = sprintf('%s/%s/%s', preproc.params.root_dir, preproc.params.output_dir, subject);
    figdir = sprintf('%s/figures', outdir);
   
    data_file = sprintf('%s/processing_results_eeg_artfrej.mat',outdir);
    
    if ~exist(data_file, 'file')
        fprintf('No data for subject %s\n', subject);
        continue;
    end
    
    load(data_file);
    
    idx_c = find(strcmp(data.eeg.ft.label,'Fz'));
    if isempty(idx_c)
        continue
    end
    
    t = data.eeg.ft.time{1};
    X = data.eeg.ft.trial{1}(idx_c,:);
    Z = (X - nanmean(X,2)) ./ nanstd(X,0,2);
    hh = plot(data.eeg.ft.time{1},Z);
    
    ax = gca;
    ax.XTick = [];
    ax.XTickLabel = [];

    ylim([-50 50]);
    xlim([min(t) max(t)]);
    
    title(subject);
    
    fprintf('Added subject %s\n', subject);
    
end

suptitle(sprintf('All time series after art rej - channel %s', channel));

end

%% Plot time series in event windows, per subject

close all;

wsize = [1400 900];
channel = 'Fz';

% Build empty subplots
clear h;
htitle = cell(8,1); hevent = cell(8,1); hselect = cell(8,1); 
hsubs = cell(8,1); i = 1;

TS = [];

h(i)=figure; h(i).Color = 'w'; h(i).Position(3:4) = wsize;
hs = zeros(length(subjects),1); 
htitle(i) = {[{'All'},{'Onset'}]};
hevent(i) = {[{'LaneChangeLeft'},{'left_change'}]};
hselect(i) = {[]};
p = 1;
for j = 1 : length(subjects)
    hs(p) = subplot(5,ceil(length(subjects)/5), p);
    p = p + 1;
end
hsubs(i) = {hs};
i = i + 1;
hs = zeros(length(subjects),1); 
h(i)=figure; h(i).Color = 'w'; h(i).Position(3:4) = wsize;
htitle(i) = {[{'All'},{'Offset'}]};
hevent(i) = {[{'LaneChangeRight'},{'right_change'}]};
hselect(i) = {[]};
p = 1;
for j = 1 : length(subjects)
    hs(p) = subplot(5,ceil(length(subjects)/5), p);
    p = p + 1;
end
hsubs(i) = {hs};
i = i + 1;
hs = zeros(length(subjects),1); 
h(i)=figure; h(i).Color = 'w'; h(i).Position(3:4) = wsize;
htitle(i) = {[{'Easy'},{'Onset'}]};
hevent(i) = {[{'LaneChangeLeft'},{'left_change'}]};
hselect(i) = {[{5},{[1 2]}]};
p = 1;
for j = 1 : length(subjects)
    hs(p) = subplot(5,ceil(length(subjects)/5), p);
    p = p + 1;
end
hsubs(i) = {hs};
i = i + 1;
hs = zeros(length(subjects),1); 
h(i)=figure; h(i).Color = 'w'; h(i).Position(3:4) = wsize;
htitle(i) = {[{'Hard'},{'Onset'}]};
hevent(i) = {[{'LaneChangeLeft'},{'left_change'}]};
hselect(i) = {[{5},{3}]};
p = 1;
for j = 1 : length(subjects)
    hs(p) = subplot(5,ceil(length(subjects)/5), p);
    p = p + 1;
end
hsubs(i) = {hs};
i = i + 1;
hs = zeros(length(subjects),1); 
h(i)=figure; h(i).Color = 'w'; h(i).Position(3:4) = wsize;
htitle(i) = {[{'Easy'},{'Offset'}]};
hevent(i) = {[{'LaneChangeRight'},{'right_change'}]};
hselect(i) = {[{5},{[1 2]}]};
p = 1;
for j = 1 : length(subjects)
    hs(p) = subplot(5,ceil(length(subjects)/5), p);
    p = p + 1;
end
hsubs(i) = {hs};
i = i + 1;
hs = zeros(length(subjects),1); 
htitle(i) = {[{'Hard'},{'Offset'}]};
hevent(i) = {[{'LaneChangeRight'},{'right_change'}]};
hselect(i) = {[{5},{3}]};
h(i)=figure; h(i).Color = 'w'; h(i).Position(3:4) = wsize;
p = 1;
for j = 1 : length(subjects)
    hs(p) = subplot(5,ceil(length(subjects)/5), p);
    p = p + 1;
end
hsubs(i) = {hs};
i = i + 1;
hs = zeros(length(subjects),1); 
h(i)=figure; h(i).Color = 'w'; h(i).Position(3:4) = wsize;
htitle(i) = {[{'Positive'},{'Onset'}]};
hevent(i) = {[{'LaneChangeLeft'},{'left_change'}]};
hselect(i) = {[{6},{[15 30]}]};
p = 1;
for j = 1 : length(subjects)
    hs(p) = subplot(5,ceil(length(subjects)/5), p);
    p = p + 1;
end
hsubs(i) = {hs};
i = i + 1;
hs = zeros(length(subjects),1); 
htitle(i) = {[{'Negative'},{'Onset'}]};
hevent(i) = {[{'LaneChangeLeft'},{'left_change'}]};
hselect(i) = {[{6},{[-25 -100]}]};
h(i)=figure; h(i).Color = 'w'; h(i).Position(3:4) = wsize;
p = 1;
for j = 1 : length(subjects)
    hs(p) = subplot(5,ceil(length(subjects)/5), p);
    p = p + 1;
end
hsubs(i) = {hs};
i = i + 1;
hs = zeros(length(subjects),1);
htitle(i) = {[{'Positive'},{'Offset'}]};
hevent(i) = {[{'LaneChangeRight'},{'right_change'}]};
hselect(i) = {[{6},{[15 30]}]};
h(i)=figure; h(i).Color = 'w'; h(i).Position(3:4) = wsize;
p = 1;
for j = 1 : length(subjects)
    hs(p) = subplot(5,ceil(length(subjects)/5), p);
    p = p + 1;
end
hsubs(i) = {hs};
i = i + 1;
hs = zeros(length(subjects),1); 
h(i)=figure; h(i).Color = 'w'; h(i).Position(3:4) = wsize;
htitle(i) = {[{'Negative'},{'Offset'}]};
hevent(i) = {[{'LaneChangeRight'},{'right_change'}]};
hselect(i) = {[{6},{[-25 -100]}]};
p = 1;
for j = 1 : length(subjects)
    hs(p) = subplot(5,ceil(length(subjects)/5), p);
    p = p + 1;
end
hsubs(i) = {hs};

xwindow = [-1000 1000];
ywindow = [-6 6];
baseline_window = [-300 0];
p = 1;
ts_started = false;

for i = 1 : length(subjects)
    
    subject = subjects{i};
    
    outdir = sprintf('%s/%s/%s', preproc.params.root_dir, preproc.params.output_dir, subject);
    figdir = sprintf('%s/figures', outdir);
    
    data_file = sprintf('%s/processing_results_eeg_artfrej.mat',outdir);
    
    if ~exist(data_file, 'file')
        fprintf('No data for subject %s\n', subject);
        p = p + 1;
        continue;
    end
    
    load(data_file);
    
    % Filter data
    cfg = [];
    cfg.hpfilter = 'yes';
    cfg.hpfreq = 1;
    cfg.hpfiltord = 2;
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 30;
    cfg.lpfiltord = 10;
    fprintf('Applying filter [%1.1f to %1.1f Hz]...', cfg.hpfreq, cfg.lpfreq);
    [~,data.eeg.ft] = evalc('ft_preprocessing(cfg, data.eeg.ft)');
    fprintf('done.\n');
    
    idx_c = find(strcmp(data.eeg.ft.label,channel));

%     if ~isempty(idx_c)
    sim = get_simulation_events_eeg( outdir, preproc, proc, data );
    trials = get_trials_eeg( data, sim, params.eeg.erp );
    Ts = 1000 / data.eeg.ft.fsample;
%     end

    % All left change trials
    for f = 1 : length(hsubs)
        hs = hsubs{f};
        figure(h(f));
        subplot(hs(p));
        axes(hs(p));
        
        if isempty(idx_c)
            title(sprintf('%s [no %s]', subject, channel));
            ax = gca;
            ax.XTick = [];
            ax.YTick = [];
            continue;
        end

        window_evt = params.eeg.erp.windows(hevent{f}{1},:);
        t = window_evt.From/Ts:window_evt.To/Ts;
        idx_bl = [find(t>=baseline_window(1),1,'first') find(t>=baseline_window(2),1,'first')];
        T = trials.(hevent{f}{2}).trl;
        if ~isempty(hselect{f})
            idx_sel = ismember(T(:,hselect{f}{1}),hselect{f}{2});
            T = T(idx_sel,:);
        end
        Ntrl = size(T,1);
        if Ntrl < 2
           warning('Not enough trials for %s, %s, %s!', subject, htitle{f}{1}, htitle{f}{2});
           title(sprintf('%s [n=%d]', subject, Ntrl));
           ax = gca;
           ax.XTick = [xwindow(1)/Ts 0 xwindow(2)/Ts];
           ax.XTickLabel = [{sprintf('%d',xwindow(1)/Ts)},{'0'},{sprintf('%d',xwindow(2)/Ts)}];
           ax.YTick = ywindow;
           xlim(xwindow/Ts);
           ylim(ywindow);
           continue;
        end
        X = zeros(Ntrl,length(t));
        for tt = 1 : Ntrl
            idx_trl = T(tt,1):T(tt,2);
            y = data.eeg.ft.trial{1}(idx_c,:);
            % Z-score this
            y = (y - nanmean(y,2)) ./ nanstd(y,0,2);
            x = y(idx_trl);
%             x = data.eeg.ft.trial{1}(idx_c,idx_trl);
            bl = y(idx_trl(idx_bl));
            x = x - mean(bl);
            X(tt,:) = x;

            % Plot lightly
            hh = plot(t,x);
            hh.Color = [.9 .9 .9];
            hh.LineWidth = 0.1;
            hold on;
        end

        % Plot mean time series
        M = nanmean(X);
        hh = plot(t,M);
        hh.Color = 'b';
        hh.LineWidth = 1.5;
        title(sprintf('%s [n=%d]', subject, Ntrl));

        ax = gca;
        ax.XTick = [xwindow(1)/Ts 0 xwindow(2)/Ts];
        ax.XTickLabel = [{sprintf('%d',xwindow(1)/Ts)},{'0'},{sprintf('%d',xwindow(2)/Ts)}];
        ax.YTick = ywindow;

        xlim(xwindow/Ts);
        ylim(ywindow);
        
    end
    
    % Add mean time series (for all channels) to TS 
    ts_started = false;
    for f = 1 : length(hsubs)
        window_evt = params.eeg.erp.windows(hevent{f}{1},:);
        t = window_evt.From/Ts:window_evt.To/Ts;
        idx_bl = [find(t>=baseline_window(1),1,'first') find(t>=baseline_window(2),1,'first')];
        T = trials.(hevent{f}{2}).trl;
        if ~isempty(hselect{f})
            idx_sel = ismember(T(:,hselect{f}{1}),hselect{f}{2});
            T = T(idx_sel,:);
        end
        Ntrl = size(T,1);
        if Ntrl < 2
           continue;
        end
        for jj = 1 : length(params.eeg.channels)
            chn = params.eeg.channels{jj};
            idx_c = find(strcmp(data.eeg.ft.label,chn));
            if isempty(idx_c)
                continue;
            end
            X = zeros(Ntrl,length(t));
            for tt = 1 : Ntrl
                idx_trl = T(tt,1):T(tt,2);
                y = data.eeg.ft.trial{1}(idx_c,:);
                % Z-score this
                y = (y - nanmean(y,2)) ./ nanstd(y,0,2);
                x = y(idx_trl);
                bl = y(idx_trl(idx_bl));
                x = x - mean(bl);
                X(tt,:) = x;
            end
            M = nanmean(X);
            if ~ts_started
                ts_started = true;
                TS.(htitle{f}{2}).(htitle{f}{1}).mean = nan(length(params.eeg.channels),length(subjects),length(M));
                TS.(htitle{f}{2}).(htitle{f}{1}).n_trials = zeros(length(subjects),1);
            end

            TS.(htitle{f}{2}).(htitle{f}{1}).mean(jj,i,:) = M;
            TS.(htitle{f}{2}).(htitle{f}{1}).n_trials(i) = Ntrl;
        end
        
    end
    
    fprintf('Added subject %s\n', subject);
    
    p = p + 1;
    
end

fprintf('Saving mean time series...');
channels = params.eeg.channels;
save(sprintf('%s/eeg_mean_tstrial.mat', proc.params.output_dir), 'TS', 'channels');
fprintf('done.\n');

fprintf('Saving figures as images...');
for i = 1 : length(h)
   figure(h(i));
   hh = suptitle(sprintf('%s Trials: Overtake %s [%s]', htitle{i}{1}, htitle{i}{2}, channel));
   hh.FontSize = 16;
   hh.FontWeight = 'bold';
   
   saveas(h(i), sprintf('%s/tstrials_%s_%s.png', fig_output_dir, htitle{i}{2}, htitle{i}{1}));
   
end

fprintf('done.\n');

%% Plot mean trial time series across subjects

if ~exist('TS','var')
   load(sprintf('%s/eeg_mean_tstrial.mat', proc.params.output_dir));
end

idx_c = find(strcmp(params.eeg.channels,channel));

% Plot onset and offset in separate figures
plot_names = [{'Onset'},{'Offset'}];
diff_names = [{'Easy'},{'Hard'}];
outcome_names = [{'Positive'},{'Negative'}];
event_names = [{'LeftLaneChange'},{'RightLaneChange'}];

for i = 1 : length(plot_names)
       
   pname = plot_names{i};
   X = squeeze(TS.(pname).All.mean(idx_c,:,:));
   
   % Plot all trials
   h = figure;
   h.Color = 'w';
   Nsubs = size(X,1);
   for k = 1 : Nsubs
       plot();
       hold on;
   end
   
    
   for j = 1 : length(diff_names)
       
       
   end
    
   for j = 1 : length(outcome_names)
       
       
   end
    
    
end




