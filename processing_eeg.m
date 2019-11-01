%% Runs statistical analyses and visualization on EEG/ET/behavioural data
%
% 0. Preprocess steps
%     a. Load data
%     b. Align time series
%     c. Load subject data to analyze
%
% 1. ERPs
%    a. Define trials/time-locking points based on simulation events
%    b. Assign conditions to trials
%    c. Compare ERPs

%% Init

do_averages = false;

addpath '/Users/lpzatr/OneDrive/OneDrive UoN/OneDrive - The University of Nottingham/synched/MATLAB/lib/fieldtrip-20180110'
addpath '/Users/lpzatr/OneDrive/OneDrive UoN/OneDrive - The University of Nottingham/synched/MATLAB/lib/cPCOH/functions';
ft_defaults

show_plots = false;

% Get params
load processing_eeg_params_osx.mat
preproc = load('preproc_params_hd.mat');
proc = load('processing_params.mat');

bad_channels = [];

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


%% Preprocess Steps

% Loads EEG and Simulation Log events and writes them to CSV files. Time series
% for these events will be zeroed to the SimulatorStarted event.
%

subjects = [{'0133'}];

for s = 1 : length(subjects)
    subject = subjects{s};
    
    fprintf('\n\n== Starting subject %s ==\n\n', subject);
    
    ok = true;
    
    outdir = sprintf('%s/%s/%s', preproc.params.root_dir, preproc.params.output_dir, subject);
    figdir = sprintf('%s/figures',outdir);
    flag_file = sprintf('%s/eeg-erp.done', outdir);
    
    if exist(flag_file,'file') && ~params.eeg.clobber
        warning('Output exists for subject %s. Skipping...', subject);
        ok = false;
    end

    %% Load data
    if ok

        zip_file = [];
         
        cfg = params.eeg.cfg;
        cfg.headerfile = sprintf('%s/%s/%s-eeg/%s.vhdr', params.eeg.data_dir, subject, subject, subject);
        if ~exist(cfg.headerfile, 'file')
           zip_file = sprintf('%s/%s/%s-eeg.zip', params.eeg.data_dir, subject, subject);
           unzip_dir = sprintf('%s/%s/%s-eeg', params.eeg.data_dir, subject, subject);
    %        fprintf('Zip file: %s\n', zip_file);
           if ~exist(zip_file, 'file')
               warning('Subject %s has no EEG data. Skipping...', subject);
               ok = false;
           else
               unzip(zip_file, unzip_dir);
               ok = exist(cfg.headerfile, 'file');
           end
        end
    end
    
    %% 
    if ok
%         cfg.dataset = sprintf('%s/%s/%s-eeg/%s.eeg', params.eeg.data_dir, subject, subject, subject);
%         [data.eeg.ft] = ft_preprocessing(cfg);
%         
%         % Remove bad channels
%         if ~isempty(bad_channels)
%             T = bad_channels{strcmp(bad_channels.Subject, subject),2};
%             if ~isempty(T)
%                 cfg2 = [];
%                 cfg2.channel = data.eeg.ft.label;
%                 for c = 1 : length(T)
%                     cfg2.channel(find(strcmp(cfg2.channel,T{c})))=[];
%                 end
%             data.eeg.ft = ft_selectdata(cfg2, data.eeg.ft);
%             fprintf('%s: Removed %d bad channels.\n', subject, length(T));
%             end
%         end
%         
%         eeg_channels = [];
%         eog_channels = [];
%         all_channels = data.eeg.ft.label;
%         idx_eog = [];
%         for c = 1 : length(data.eeg.ft.label)
%            channel = data.eeg.ft.label{c};
%            if contains(channel,'EOG')
%                eog_channels = [eog_channels {channel}];
%                idx_eog = [idx_eog c];
%            else
%                eeg_channels = [eeg_channels {channel}];
%            end
%         end
%         
%         marker_file = sprintf('%s/%s/%s-eeg/%s.vmrk', params.eeg.data_dir, subject, subject, subject);
%         cfg.event = ft_read_event(marker_file);
%         idx_eeg_start = find(strcmp(vertcat({cfg.event.value}),'S128'));
%         idx_eeg_start = idx_eeg_start(1);
%         t_eeg_start = data.eeg.ft.time{1}(idx_eeg_start);
%         
%         preproc = load('preproc_params_hd.mat');
% 
%         % -- First do EEG events --
% 
%         is_stim = strcmp({cfg.event(:).type}, 'Stimulus');
%         tidx = cell2mat({cfg.event(is_stim).sample});
%         value = {cfg.event(is_stim).value};
%         type = {cfg.event(is_stim).type};
%         triggers = zeros(1,length(value));
% 
%         for i = 1 : length(value)
%             event = value{i};
%             triggers(i) = str2num(event(2:end));
%         end
% 
%         % First trigger is simulation start
%         start_sim = tidx(1);
%         t_simstart_eeg = data.eeg.ft.time{1}(start_sim);
%         triggers = triggers(2:end);
%         value = value(2:end);
%         type = type(2:end);
% 
%         % Add 256 to each reset of byte values (SerialByte = 1)
%         to_add = 0;
%         first_found = false;
%         for i = 1 : length(triggers)
%             if triggers(i) == 1
%                 if first_found
%                     to_add = to_add + 256;
%                 else
%                     first_found = true;
%                 end
%             end
%             triggers(i) = triggers(i) + to_add;
%         end
% 
%         times = data.eeg.ft.time{1}(tidx);
%         times = times - t_simstart_eeg;
%         times = times(2:end);
%         tidx = tidx(2:end);
% 
%         data.eeg.t_simstart = t_simstart_eeg;
%         data.eeg.trigger_idx = tidx;
%         data.eeg.cfg = cfg;
%         data.eeg.events = table(times', tidx', triggers', value', type');
%         data.eeg.events.Properties.VariableNames = {'Time', 'Index', 'Trigger', 'Value', 'Type'};
% 
%         [t_delta, t_simstart_sim] = get_simstart_offset( preproc.params, subject ); % Unit is ms
%         data.eeg.ft.time = {data.eeg.ft.time{1} - t_simstart_eeg};

        
        data = load_eeg_data( params, subject );

         %% Perform ICA - Remove artifactual components
        
        if params.eeg.ica.apply
            
            load(sprintf('%s/processing_results_eeg_ica.mat',outdir), 'data_flt');
            
%             cfg_ica = data.eeg.cfg;
%             cfg_ica.continuous = 'yes';
%             cfg_ica.hpfilter = 'yes';
%             cfg_ica.hpfreq = 1;
%             cfg_ica.hpfiltord = 6;
%             cfg_ica.lpfilter = 'yes';
%             cfg_ica.lpfreq = 60;
%             cfg_ica.lpfiltord = 6;
%             cfg_ica.demean = 'no';
%             cfg_ica.reref = 'yes';
%             cfg_ica.refchannel = 'all';
%             
%             data_ica = ft_preprocessing(cfg_ica);
%             dt_eeg = data.eeg.events{1,1};
%             data_ica.time = {data.eeg.ft.time{1} - dt_eeg};
%             
% %             cfg_ica = [];
% %             cfg_ica.resamplefs = 300;
% %             cfg_ica.detrend    = 'no';
% %             data_ica = ft_resampledata(cfg_ica, data_ica);
%             
%             cfg_ica = [];
%             cfg_ica.method = 'runica'; % this is the default and uses the implementation from EEGLAB
% 
%             results.eeg.ica.comp = ft_componentanalysis(cfg_ica, data_ica);
%             
%             cfg_ica = [];
%             cfg_ica.colormap = 'jet';
%             cfg_ica.component = 1:20;       % specify the component(s) that should be plotted
%             cfg_ica.layout    = 'acticap-64ch-standard2.mat'; % specify the layout file that should be used for plotting
%             cfg_ica.comment   = 'no';
%             h = figure;
%             h.Color = 'w';
%             ft_topoplotIC(cfg_ica, results.eeg.ica.comp);
%             resize_window(h,[1000 1000]);
%             saveas(h,sprintf('%s/ica_topoplot.png',figdir));
% 
%             cfg_ica.viewmode = 'component';
%             ft_databrowser(cfg_ica, results.eeg.ica.comp);
%             resize_window(gcf,[1500 900]);
%             saveas(gcf,sprintf('%s/ica_browser.fig',figdir));
% 
%             % Prompt for bad channels
%             channels = input('Enter ICA components to remove: ','s');
%             to_rem = cellfun(@str2num, strsplit(channels, ' '));
%             
%             if length(to_rem) > 0
%                fprintf('Removing %d components...\n', length(to_rem));
%                cfg_ica = [];
%                cfg_ica.component = to_rem;
%                data_flt = ft_rejectcomponent(cfg_ica, results.eeg.ica.comp, data_ica);
%                results.eeg.ica.removed = to_rem;
%             end
%             
%             fprintf('Done. Saving results...\n');
%             save(sprintf('%s/processing_results_eeg_ica.mat',outdir), 'results', 'data_flt');
%             
%             close all;

        end
        
        %% Define simulation events on this timeline

        %Load lange change events
        input_file = sprintf('%s/events-LaneChange.csv', outdir);
        [values, hdr] = import_log(input_file, preproc.params.simlog.lanechange_format);

        T = values{1};
        mylog.times = double(T - t_simstart_sim) / 1000; % Unit is s
        mylog.triggers = values{find(strcmp(hdr,'AdjSerialByte'))};
        mylog.logid = values{find(strcmp(hdr,'LogId'))};
        mylog.type = values{find(strcmp(hdr,'EventType'))};
        mylog.lane_from = values{find(strcmp(hdr,'LaneFrom'))};
        mylog.lane_to = values{find(strcmp(hdr,'LaneTo'))};
        mylog.sim_time =values{find(strcmp(hdr,'SimulationTime'))};
        mylog.cycle=values{find(strcmp(hdr,'Sim:Game:Cycle'))};
        mylog.repeat=values{find(strcmp(hdr,'Sim:Game:Repeat'))};

        sim.events = table(mylog.times, mylog.triggers, mylog.logid, mylog.type, mylog.lane_from, ...
                           mylog.lane_to, mylog.sim_time, mylog.cycle, mylog.repeat);
        sim.events.Properties.VariableNames = {'Time', 'Trigger', 'LogId', 'Type', 'LaneFrom', ...
                                               'LaneTo', 'SimTime', 'Cycle', 'Repeat'};
                              
        writetable(sim.events, sprintf('%s/sim-events-lanechange.csv', outdir));
        
        %% Define passing and baseline epochs
        proc_results_file = sprintf('%s/processing_results.mat', outdir);
        proc_results = load(proc_results_file);
        
        sim.epochs.idx_baseline = proc_results.results.epochs.intervals.baseline.idx;
        sim.epochs.idx_passing = proc_results.results.epochs.intervals.passing.idx;
        
        min_epoch = 500;
        
        trl_epochs = zeros(0,5);
        for b = 1 : size(sim.epochs.idx_baseline,1)
           idxb =  sim.epochs.idx_baseline(b,:);
           if idxb(2) - idxb(1) >= min_epoch
                trl_epochs(end+1,:) = [idxb(1),idxb(2),0,1,1];
           end
        end
        for b = 1 : size(sim.epochs.idx_passing,1)
           idxb =  sim.epochs.idx_passing(b,:);
           if idxb(2) - idxb(1) >= min_epoch
                trl_epochs(end+1,:) = [idxb(1),idxb(2),0,1,2];
           end
        end
        

        %% Load eye data 
        et_results_file = sprintf('%s/%s/%s/results.mat',preproc.params.root_dir, ...
                                  preproc.params.output_dir, subject);
        et_results = load(et_results_file); 
        
        dt_eye = single(et_results.data.eye.log.messages{1}(1)/1000 - et_results.data.eye.t_start);
        % Align data to first log events
        et_results.results.t = et_results.results.t - dt_eye;
        dt_eeg = data.eeg.events{1,1};
        
        cfg = data.eeg.cfg;
        cfg.continuous = 'yes';
        cfg.hpfilter = 'yes';
        cfg.hpfreq = 0.5;
        cfg.hpfiltord = 6;
        cfg.lpfilter = 'yes';
        cfg.lpfreq = 60;
        cfg.lpfiltord = 6;
        cfg.demean = 'no';
%         cfg.reref = 'yes';
%         cfg.refchannel = 'all';
        
        % If ICA wasn't run, load raw EEG
        if ~exist('data_flt','var')
            data_flt = ft_preprocessing(cfg);
            data_flt.time = {data.eeg.ft.time{1} - dt_eeg};
        end
        
        % Artifact rejection based on ET-defined saccades and eyeblinks
        if ~strcmp(params.eeg.artifacts.reject,'none')
           t_eye = et_results.results.t / 1000;
           t_eeg = data_flt.time{1};
           
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
           for i = 1 : size(idx_sacc,1)
                xs = t_eye(idx_sacc(i,:));
                idx0 = findInSorted(t_eeg, xs(1));
                idx1 = findInSorted(t_eeg, xs(2));
                artifacts = [artifacts;[idx0 idx1]];
           end
           
           cfg_art = [];
           cfg_art.continuous = 'yes';
           cfg_art.artfctdef.zvalue.channel = eeg_channels;
           cfg_art.artfctdef.zvalue.cutoff = 20;
           cfg_art.artfctdef.zvalue.trlpadding = 0;
           cfg_art.artfctdef.zvalue.fltpadding = 0;
           cfg_art.artfctdef.zvalue.artpadding = 0.05;
           [~, artifactz] = ft_artifact_zvalue(cfg_art, data_flt);
           
           cfg_art = [];
           cfg_art.artfctdef.reject          = 'nan'; %params.eeg.artifacts.reject;
           cfg_art.artfctdef.eog.artifact    = artifacts;
           cfg_art.artfctdef.zvalue.artifact = artifactz;
           cfg_art.artfctdef.minaccepttim    = params.eeg.artifacts.minaccepttim_fft;
           data.eeg.cfg.artfctdef = cfg_art.artfctdef;
           data_flt = ft_rejectartifact(cfg_art, data_flt);
           
        end
        
        save(sprintf('%s/processing_results_eeg_artfrej.mat',outdir), 'data_flt', '-v7.3');
        
        %if show_plots
           plot_eeg_blinks(et_results.results, data_flt, [{'vEOGover'},{'Fz'},{'FC5'},{'POz'},{'T7'}], 6, ...
               sprintf('%s/eeg_on_blinks.png', figdir) );
        %end
        
%         clear data_flt;


        %% Define trials from events
        %
        % Note in Fieldtrip, trials are defined by indices

        W = params.eeg.erp.windows;

        % LaneChangeLeft events
        T = sim.events(strcmp(sim.events.LaneFrom,'Lane.hiway.1'),{'Time','Trigger'});
        window = W('LaneChangeLeft',:);
        trial_idx = window.TrialIdx;

        % Overtake difficulty
        process_results_file = sprintf('%s/%s/processing_results.mat',proc.params.data_dir,subject);
        processing_results = load(process_results_file);
        diffs = processing_results.results.events.left_change.diffs;

        Ts = 1 / data.eeg.ft.fsample; % Sample period in ms
        window = [window.From / Ts, window.To / Ts];
        trl = zeros(0,5);
        trl_baseline = zeros(0,5);
        baseline_offset = -20.0; % seconds
        baseline_offset = round(baseline_offset * 1000 / Ts);

        time_tol = 0.01;

        % Define Left-lane-change trials
        for i = 1 : height(T)
            record_i = data.eeg.events(data.eeg.events.Trigger==T.Trigger(i),:);
            if isempty(record_i)
                warning('   No EEG event matching trigger %d...', T.Trigger(i));
            else
                ok=1;
                % Deal with double triggers (should not occur but do)
                if height(record_i) > 1
                    t_eeg = record_i{:,{'Time'}};
                    t_sim = T{i,{'Time'}};
                    ok = abs(t_eeg - t_sim) < time_tol;
                end
                for j = 1 : length(ok)
                    if ok(j)
                        record_ij = record_i(j,:);
                        window_ij = [record_ij.Index + window(1), record_ij.Index + window(2)];
                        trl(end+1,:) = [window_ij(1), window_ij(2), window(1), trial_idx, diffs(i)]; 
                        bw_window = window_ij + baseline_offset;
                        if bw_window(1) > 0 && bw_window(2) < length(data.eeg.ft.time{1})
                            trl_baseline(end+1,:) = [bw_window(1), bw_window(2), window(1), trial_idx, diffs(i)];
                        else
                            trl_baseline(end+1,:) = nan(5,1);
                        end
                    end
                end
            end
        end
        
        idx_nan = find(any(isnan(trl_baseline),2));
        if ~isempty(idx_nan)
           trl(idx_nan,:) = [];
           trl_baseline(idx_nan,:) = [];
        end
        
        
        %% Hilbert transforms
        
        if params.eeg.hilbert.apply
        
            cfg_hb = data.eeg.cfg;
            cfg_hb.continuous = 1;
            cfg_hb.demean = 'no';

            T = params.eeg.hilbert.bands;
            
            N_bands = height(T);
            results.eeg.hilbert.bands = T;

            L = 0;
            N_chan = length(data_flt.label);
            for j = 1 : length(data_flt.trial)
                L = L + length(data_flt.time{j});
            end

            cfg_interp = [];
            cfg_interp.method = 'linear';
            cfg_interp.prewindow = 0.004;
            cfg_interp.postwindow = 0.004;
            cfg_interp.feedback = 'no';

            idata = ft_interpolatenan_corrected(cfg_interp, data_flt);
            
            for bb = 1 : N_bands
                
                freqband = [T.From(bb) T.To(bb)];
                fprintf('\nComputing envelopes for %s: [%1.1f to %1.1f Hz]\n', ...
                            T.Band{bb}, freqband(1), freqband(2));

                % make filter
                cfg_hb = [];
                cfg_hb.filtord = T.Filtord(bb);
                tol = 100;
                [b,a] = butter(cfg_hb.filtord,2*freqband/idata.fsample, 'bandpass');

                h = figure('visible','off');
%                 fvtool(b,a, 'Fs', idata.fsample);
                fmin = max(0.1,freqband(1)-10);
                fmax = max(20, freqband(2)+10);
                freqz(b,a,fmin:(fmax-fmin)/1000:fmax,idata.fsample);
                saveas(h,sprintf('%s/hilbert_filter_%s.png', figdir, T.Band{bb}));
                close(h);
                
                % apply filter to data
                N_t = size(idata.trial{1},2);
                datfilt = zeros(length(eeg_channels),N_t);
                envelopes = zeros(length(eeg_channels),N_t);
                for cc = 1 : length(eeg_channels)
                    idx = find(strcmp(idata.label, eeg_channels{cc}));
                    X = idata.trial{1}(idx,:);
                    Y = filtfilt(b,a,X);
                    datfilt(cc,:) = Y;
                    envelopes(cc,:) = hilbert(Y);
                end
                
                results.eeg.hilbert.cfg(bb) = cfg_hb;
                results.eeg.hilbert.Fs = idata.fsample;
                results.eeg.hilbert.filtered(bb) = {datfilt};
                results.eeg.hilbert.channels = eeg_channels;
                results.eeg.hilbert.envelopes(bb) = {envelopes};
                
            end

           plot_params = preproc.params.plots.events;
           plot_params.plot_overtakes = true;
           plot_params.plot_saccades = true;
           plot_params.xlim=[0 0.5];

           plot_params.eeg.stdev = 6;
           plot_params.eeg.line_widths = [0.5 1.0 2];
           plot_params.eeg.alpha = [0.3, 0.6, 1.0];
           plot_params.eeg.height = 20;
           plot_params.eeg.interpolate = true;
           plot_params.eeg.smooth = 0;

           for bb = 1 : N_bands

               plot_params.to_file =  [{sprintf('%s/hilbert_ts_%s.fig', figdir, T.Band{bb})} ...
                                       {sprintf('%s/hilbert_ts_%s.png', figdir, T.Band{bb})}];

               data_f = data_flt;
               data_f.trial={abs(results.eeg.hilbert.filtered{bb})};
               data_f.label = results.eeg.hilbert.channels;
               data_h = data_flt;
               data_h.trial={abs(results.eeg.hilbert.envelopes{bb})};
               data_h.label = results.eeg.hilbert.channels;

               plot_events2(et_results.results, plot_params, [{'Rounds'},{'LaneChanges'},{'Baseline'},{'Saccades'},{'Overtakes'}], ...
                            [{data_flt}, {data_f}, {data_h}], [{'AF7'},{'Fz'},{'FC5'},{'Cz'},{'T7'}]);
               hh = title(sprintf('Hilbert envelopes for %s [%1.1f to %1.1f Hz]', T.Band{bb}, T.From(bb), T.To(bb)));
               hh.FontSize = 16;
           end

            save(sprintf('%s/processing_results_eeg_hilbert.mat',outdir), 'results', '-v7.3');
            clear data_hilbert data_f data_h;
        
        end
        
        %% Run analyses

        cfg = data.eeg.cfg;
        cfg.continuous = 0;
        cfg.hpfilter = 'yes';
        cfg.hpfreq = 0.5;
        cfg.hpfiltord = 2;
        cfg.lpfilter = 'yes';
        cfg.lpfreq = 60;
        cfg.lpfiltord = 10;
        cfg.demean = 'no';
        cfg.baselinewindow = [-0.3 -0.1];
        cfg.pad='nextpow2';
                
        % For resampling
        cfg2 = [];
        cfg2.resamplefs = 256;
        cfg2.detrend = 'no';
        
        % Apply trials to the data
        cfg.trl = trl;
        cfg = ft_rejectartifact(cfg);
        data_erp.left_change.all = ft_preprocessing(cfg, data_flt);
        data_erp.left_change.all = ft_resampledata(cfg2, data_erp.left_change.all);

        fprintf('\n == Loading data: Easy v. Difficult ==\n\n');
        
        % Easy
        idx = find(trl(:,5)==1);
        cfg.trl = trl(idx,1:4);
        cfg = ft_rejectartifact(cfg, data_flt);
        data_erp.left_change.easy = ft_preprocessing(cfg, data_flt);
        data_erp.left_change.easy = ft_resampledata(cfg2, data_erp.left_change.easy);

        % Difficult
        idx = find(trl(:,5)==2);
        cfg.trl = trl(idx,1:4);
        cfg = ft_rejectartifact(cfg, data_flt);
        data_erp.left_change.difficult = ft_preprocessing(cfg, data_flt);
        data_erp.left_change.difficult = ft_resampledata(cfg2, data_erp.left_change.difficult);

        fprintf('\n == Loading data: Baseline v. Passing ==\n\n');
        
        % Baseline
        idx = find(trl_epochs(:,5)==1);
        cfg.trl = trl_epochs(idx,1:4);
        cfg = ft_rejectartifact(cfg, data_flt);
        data_epochs.baseline = ft_preprocessing(cfg, data_flt);
        data_epochs.baseline = ft_resampledata(cfg2, data_epochs.baseline);
        
        % Passing
        idx = find(trl_epochs(:,5)==2);
        cfg.trl = trl_epochs(idx,1:4);
        cfg = ft_rejectartifact(cfg, data_flt);
        data_epochs.passing = ft_preprocessing(cfg, data_flt);
        data_epochs.passing = ft_resampledata(cfg2, data_epochs.passing);

        
        %% ERP Analysis

        if params.eeg.erp.apply
        
        fprintf('\n == ERP: Passing events ==\n\n');
        
            cfg_erp = [];
            cfg_erp.trials = 'all';
            cfg_erp.vartrllength = 2;

            data_erp.left_change.all = ft_timelockanalysis(cfg_erp, data_erp.left_change.all);
            data_erp.left_change.easy = ft_timelockanalysis(cfg_erp, data_erp.left_change.easy);
            data_erp.left_change.difficult = ft_timelockanalysis(cfg_erp, data_erp.left_change.difficult);

            save(sprintf('%s/processing_results_eeg_erp.mat',outdir), 'data_erp', '-v7.3');
            clear data_erp;
        
        end
        
        %% FFT & Time/Freq Analyses - Baseline versus Passing Epochs
        
        fprintf('\n == FFT: Baseline v. Passing ==\n\n');
        
        cfg_fft = [];
        cfg_fft.method     = 'mtmfft';
        cfg_fft.output     = 'pow';
        cfg_fft.channel    = eeg_channels;
        cfg_fft.keeptrials = 'yes';
        cfg_fft.taper      = 'hanning';
        cfg_fft.foilim     = [cfg.hpfreq cfg.lpfreq];
%         cfg_fft.foi     = cfg.hpfreq : 0.1 : cfg.lpfreq;
        
        % Baseline
        data_fft.epochs.baseline = ft_freqanalysis(cfg_fft, data_epochs.baseline);

        % Passing
        data_fft.epochs.passing = ft_freqanalysis(cfg_fft, data_epochs.passing);
        
        save(sprintf('%s/processing_results_eeg_epochs.mat',outdir), 'data_epochs', '-v7.3');
%         clear data_epochs;
               
        %% FFT
        
        if params.eeg.fft.apply
        
            cfg_fft = params.eeg.fft.cfg;
%             cfg_fft.method     = 'mtmfft';
%             cfg_fft.output     = 'pow';
%             cfg_fft.channel    = eeg_channels;
%             cfg_fft.keeptrials = 'yes';
%             cfg_fft.taper      = 'hanning';
%             cfg_fft.pad =      'nextpow2';

            fstep = 1 / ((window(2) - window(1)) / 1000);
            tstep = Ts * 10;

            cfg_timefreq = [];
            cfg_timefreq.method     = 'wavelet';
            cfg_timefreq.output     = 'fourier';
            cfg_timefreq.taper      = 'hanning';
            cfg_timefreq.channel    = eeg_channels;
            cfg_timefreq.keeptrials = 'yes';
            cfg_timefreq.keeptapers = 'yes';
            cfg_timefreq.foi        = 1:fstep:40;
    %         cfg_timefreq.t_ftimwin  = ones(length(cfg_timefreq.foi),1).*0.5;
            cfg_timefreq.toi        = [window(1):tstep:window(2)] / 1000;
            cfg_timefreq.pad =      'nextpow2';

            % FFT: All
            fprintf('\n == FFT: Easy v. Difficult ==\n\n');
            cfg.artfctdef.minaccepttim   = params.eeg.artifacts.minaccepttim_fft;

            cfg.trl = trl(:,1:4);
            mydata = ft_preprocessing(cfg);
            mydata = ft_rejectartifact(cfg,  mydata);
            data_fft.left_change.all = ft_freqanalysis(cfg_fft, mydata);

            % FFT: Easy
            idx = find(trl(:,5)==1);
            cfg.trl = trl(idx,1:4);
            mydata = ft_preprocessing(cfg);
            mydata = ft_rejectartifact(cfg,  mydata);
            data_fft.left_change.easy = ft_freqanalysis(cfg_fft, mydata);

            % FFT: Difficult
            idx = find(trl(:,5)==2);
            cfg.trl = trl(idx,1:4);
            mydata = ft_preprocessing(cfg);
            mydata = ft_rejectartifact(cfg,  mydata);
            data_fft.left_change.difficult = ft_freqanalysis(cfg_fft, mydata);

            save(sprintf('%s/processing_results_eeg_fft.mat',outdir), 'data_fft', '-v7.3');
            clear data_fft;
        
        end
        
        %% Time/Freq Analyses - Easy versus Difficult Events
        
        if params.eeg.timefreq.apply
        
            % Time-freq: All
            fprintf('\n == Time-freq: Easy v. Difficult ==\n\n');
            cfg = data.eeg.cfg;
            cfg.continuous = 0;
            cfg.hpfilter = 'yes';
            cfg.hpfreq = 0.5;
            cfg.hpfiltord = 2;
            cfg.lpfilter = 'yes';
            cfg.lpfreq = 60;
            cfg.lpfiltord = 10;
            cfg.demean = 'no';
            cfg.baselinewindow = [-0.3 -0.1];
            cfg.pad='nextpow2';
        
            cfg.artfctdef.eog.artifact   = artifacts;
            cfg.artfctdef.minaccepttim   = params.eeg.artifacts.minaccepttim_timefreq;
            cfg.artfctdef.reject = 'nan';
            cfg.detrend = 'yes';

            cfg_interp = [];
            cfg_interp.method = 'linear';
            cfg_interp.prewindow = 0.004;
            cfg_interp.postwindow = 0.004;
            cfg_interp.feedback = 'no';

            cfg.trl = trl(:,1:4);
            mydata = ft_preprocessing(cfg);
            mydata = ft_rejectartifact(cfg,  mydata);
            mydata = ft_interpolatenan_corrected(cfg_interp, mydata);

            cfg_bl = cfg;
            cfg_bl.trl = trl_baseline(:,1:4);
            bldata = ft_preprocessing(cfg_bl); 
            bldata = ft_rejectartifact(cfg_bl,  bldata);
            bldata = ft_interpolatenan_corrected(cfg_interp, bldata);

            mydata.event = cell(1,size(trl,1));
            bldata.event = cell(1,size(trl,1));
            for a = 1 : size(trl,1)
                mydata.event(a) = {0};
                bldata.event(a) = {0};
            end

            [~,idx_eeg] = intersect(data_flt.label, eeg_channels);

            cfg_timefreq = params.eeg.timefreq.cfg;
%             cfg_timefreq.foi = [2:0.5:14,15:24,26:2:40]; % Hz
%             cfg_timefreq.foi = [2:0.5:14,15:22,24:2:30];
%             cfg_timefreq.wavelet = 'cmor6-1';
%             cfg_timefreq.dt = 0.002; %s
%             cfg_timefreq.toi = [-15,10]; %s
%             cfg_timefreq.Zscore = true;

            cfg_timefreq.channels = idx_eeg;
            
    %         cfs_bl = getSpectra_baseline(cfg_timefreq, bldata, []);
            
            N_trial = length(mydata.trial);
    
            if cfg_timefreq.heavy
                cfs = cell(length(idx_eeg),1);
                for c = 1 : length(idx_eeg)
                    idx_c = idx_eeg(c);
                    fprintf(' Computing spectra for channel %s...', mydata.label{idx_c});
                    mydata_c = mydata;
                    bldata_c = bldata;
                    mydata_c.label = mydata.label(idx_c);
                    bldata_c.label = bldata_c.label(idx_c);
                    for jj = 1 : N_trial
                        mydata_c.trial(jj) = {mydata_c.trial{jj}(idx_c,:)};
                        bldata_c.trial(jj) = {bldata_c.trial{jj}(idx_c,:)};
                    end
                    cfg_c = cfg_timefreq;
                    cfg_c.channels = 1;
                    cfs_c = getSpectra_baseline(cfg_c, mydata_c, bldata_c);
                    cfs(c) = {cfs_c};
                    fprintf('done.\n');
                end
            else
                cfs = getSpectra_baseline(cfg_timefreq, mydata, bldata);
            end
    %        
            data_timefreq = [];
            data_timefreq.left_change = cell(length(idx_eeg),1);

            for c = 1 : length(idx_eeg)
                % All compared to baseline
                if cfg_timefreq.heavy
                    X_all = squeeze(abs(cfs{c}.powspctrm(:,:,1,:,:)));
                else
                    X_all = squeeze(abs(cfs.powspctrm(:,:,c,:,:)));
                end
                s_all = squeeze(nanstd(X_all,1,1));
                X_all = squeeze(nanmean(X_all,1));
                n_all = size(trl,1);
  
                df = n_all-1;
                s_p = s_all; %sqrt((s_all.^2 + s_bl.^2)/2);
                C_all = X_all;
                T_all = C_all./(s_p * sqrt(2/n_all));
                P_all = 2 * (1 - tcdf(abs(T_all), df));
                
                % Time-freq: Easy
                idx = find(trl(:,5)==1);
                n_easy = length(idx);
%                 X_easy = abs(squeeze(cfs.powspctrm(:,idx,c,:,:)));
                if cfg_timefreq.heavy
                    X_easy = squeeze(abs(cfs{c}.powspctrm(:,idx,1,:,:)));
                else
                    X_easy = squeeze(abs(cfs.powspctrm(:,idx,c,:,:)));
                end
                s_easy = squeeze(nanstd(X_easy,1,1));
                X_easy = squeeze(nanmean(X_easy,1));

        %         % Time-freq: Difficult
                idx = find(trl(:,5)==2);
                n_diff = length(idx);
%                 X_diff = abs(squeeze(cfs.powspctrm(:,idx,c,:,:)));
                if cfg_timefreq.heavy
                    X_diff = squeeze(abs(cfs{c}.powspctrm(:,idx,1,:,:)));
                else
                    X_diff = squeeze(abs(cfs.powspctrm(:,idx,c,:,:)));
                end
                s_diff = squeeze(nanstd(X_diff,1,1));
                X_diff = squeeze(nanmean(X_diff,1));

                df = n_easy+n_diff-2;
                s_p = sqrt(((n_easy-1)*s_easy.^2 + (n_diff-1)*s_diff.^2)/df);
                C_de = X_diff - X_easy;
                T_de = C_de./(s_p * sqrt(1/n_easy + 1/n_diff));
                P_de = 2 * (1 - tcdf(abs(T_de), df)); 

                ch_result = [];
                ch_result.cfg = cfg_timefreq;
                ch_result.channels = data_flt.label(idx_eeg);
                ch_result.stats.all_baseline.c = C_all;
                ch_result.stats.all_baseline.t = T_all;
                ch_result.stats.all_baseline.p = P_all;
                ch_result.stats.all_baseline.n_all = n_all;
                ch_result.stats.easy_diff.c = C_de;
                ch_result.stats.easy_diff.t = T_de;
                ch_result.stats.easy_diff.p = P_de;
                ch_result.stats.easy_diff.n_diff = n_diff;
                ch_result.stats.easy_diff.n_easy = n_easy;

                data_timefreq.left_change(c) = {ch_result};

                clear mydata bldata;
            end

            save(sprintf('%s/processing_results_eeg_timefreq.mat',outdir), 'data_timefreq', '-v7.3');
            clear data_timefreq;
        
        end
        
        %% Save and Tidy
        remfile = sprintf('%s/processing_results_eeg.mat',outdir);
        if exist(sprintf('%s/processing_results_eeg.mat',outdir), 'file')
           delete(remfile);
        end    
        
        % Clean up by removing unzipped data
        if ~isempty(zip_file) && exist(zip_file, 'file')
            if exist(unzip_dir, 'dir')
                rmdir(unzip_dir, 's');
%                 fprintf('  *DEBUG: Removed unzip data for %s\n', subject);
            end
        end

        clear data;
        clear cfg;
        
        fid = fopen(flag_file, 'w');
        fclose(fid);
        
        fprintf('Done processing subject %s\n', subject);
    end
    
end


%% Average across subjects

if do_averages

    clear data;

    load(proc.params.qc.file);
    qc_score = cell2mat(qc_eeg(:,1));
    idx = qc_score>1; %=proc.params.qc.cutoff;
    subjects = qc_eeg(idx,2);

    tldata.all = [];
    tldata.easy = [];
    tldata.difficult = [];

    fftdata.all = [];
    fftdata.easy = [];
    fftdata.difficult = [];

    epochdata.baseline = [];
    epochdata.passing = [];

    for s = 1 : length(subjects)
        subject = subjects{s};
        ok = true;

        outdir = sprintf('%s/%s/%s', preproc.params.root_dir, preproc.params.output_dir, subject);
        flag_file = sprintf('%s/eeg-erp.done', outdir);

        if ~exist(flag_file,'file')
            warning('Skipping subject %s...', subject);
            ok = false;
        end

        if ok
            results_file = sprintf('%s/processing_results_eeg_erp.mat',outdir);
            load(results_file);

            tldata.all = [tldata.all {data_erp.left_change.all}];
            tldata.easy = [tldata.easy {data_erp.left_change.easy}];
            tldata.difficult = [tldata.difficult {data_erp.left_change.difficult}];
            clear data_erp;

            results_file = sprintf('%s/processing_results_eeg_timefreq.mat',outdir);
            load(results_file);
    %         fftdata.all = [fftdata.all {data_timefreq.left_change.all}];
            fftdata.easy = [fftdata.easy {data_timefreq.left_change.easy}];
            fftdata.difficult = [fftdata.difficult {data_timefreq.left_change.difficult}];
            clear data_timefreq;

            results_file = sprintf('%s/processing_results_eeg_epochs.mat',outdir);
            load(results_file);
            epochdata.baseline = [epochdata.baseline {data_epochs.baseline}];
            epochdata.passing = [epochdata.passing {data_epochs.passing}];
            clear data_epochs;

            fprintf('Added subject %s.\n', subject);
        end

    end

    clear data;

    if ~exist(params.eeg.output_dir, 'dir')
       mkdir(params.eeg.output_dir);
    end

    save(sprintf('%s/erp-left_change-difficulty.mat',params.eeg.output_dir), 'tldata');
    save(sprintf('%s/tlock-left_change-difficulty.mat',params.eeg.output_dir), 'fftdata', '-v7.3');
    save(sprintf('%s/fft-baseline-passing.mat',params.eeg.output_dir), 'epochdata', '-v7.3');


    %% Grand averages

    if ~exist('tldata','var')
        load(sprintf('%s/erp-left_change-difficulty.mat',params.eeg.output_dir));
    end

    cfg = [];
    cfg.channel   = 'all';
    cfg.latency   = 'all';
    cfg.parameter = 'avg';

    gavr.all = ft_timelockgrandaverage(cfg, tldata.all{:});
    gavr.easy = ft_timelockgrandaverage(cfg, tldata.easy{:});
    gavr.difficult = ft_timelockgrandaverage(cfg, tldata.difficult{:});

    if ~exist('fftdata','var')
        load(sprintf('%s/tlock-left_change-difficulty.mat',params.eeg.output_dir));
    end

    clear fft_gavr;
    X = [fftdata.all{:}];
    X_easy = [fftdata.easy{:}];
    X_diff = [fftdata.difficult{:}];

    N_subj = length(subjects);
    N_frq = length(X(1).freq);
    N_t = length(X(1).time);
    N_ch = length(params.eeg.channels);

    Y_all = nan(N_subj, N_ch, N_frq, N_t);
    Y_easy = nan(N_subj, N_ch, N_frq, N_t);
    Y_diff = nan(N_subj, N_ch, N_frq, N_t);

    for i = 1 : N_subj

        chn = X(i).label;
        idx_c = [];
        for j = 1 : length(chn)
            if find(strcmp(params.eeg.channels, chn(j)))
                idx_c = [idx_c find(strcmp(params.eeg.channels,chn(j)))]; 
            else
                fprintf('%s\n', chn{j});
            end
        end

        Y_all(i,idx_c,:,:) = X(i).powspctrm;
        Y_easy(i,idx_c,:,:) = X_easy(i).powspctrm;
        Y_diff(i,idx_c,:,:) = X_diff(i).powspctrm;

    end



    %% Plot electrodes ERP

    figure;
    cfg = [];
    cfg.layout = 'acticap-64ch-standard2.mat';
    cfg.interactive = 'yes';
    cfg.showoutline = 'yes';
    cfg.ylim=[-4 4];
    cfg.xlim=[-0.5 1.5];
    cfg.showlabels='yes';
    cfg.box='yes';

    h=ft_multiplotER(cfg, gavr.easy, gavr.difficult);


    %% Single channel plots 

    show_channels = [{'Fz'},{'F8'},{'AF7'}];

    for i = 1 : length(show_channels)

        chan = show_channels{i};
        idx = find(strcmp(gavr.easy.cfg.channel, chan));

        avrvals = zeros(size(gavr.easy.avg,2),2);
        avrvals(:,1) = gavr.easy.avg(idx,:);
        avrvals(:,2) = gavr.difficult.avg(idx,:);

        t = 1000 * repmat(gavr.easy.time, [2 1])';

        h = figure;
        set(h, 'Color', 'w');
        hh = plot(t, avrvals);

        set(hh,'LineWidth', 1.3);
        hh = title(sprintf('ERPs time-locked to left lane change [%s]', chan));
        set(hh,'FontSize', 16);
        hh = xlabel('Time rel. to button press (ms)');
        set(hh,'FontSize', 14);
        hh = ylabel('Potential (mV)');
        set(hh,'FontSize', 14);

        ylim([-4 4]);

        hold on;
        hh = plot([0 0],[-4 4],'k');
        set(hh,'LineWidth',1.2);

        legend([{'Easy'},{'Difficult'}]);

    end


    %% Plot electrodes FFT



    %% Plot single channels - FFT

    freq = fftdata.all{1}.freq;
    idx_freq = find(freq >= 2 & freq <= 50);

    show_channels = [{'Fz'},{'F8'},{'AF7'}];
    zstar = 1.96;

    for i = 1 : length(show_channels)

        chan = show_channels{i};
        idx = find(strcmp(eeg_channels, chan));
        avrvals = zeros(size(fft_gavr.easy.mean,2),2);
        civals = zeros(size(fft_gavr.easy.mean,2),2);
        avrvals(:,1) = fft_gavr.easy.mean(idx,:);
        avrvals(:,2) = fft_gavr.difficult.mean(idx,:);
        civals(:,1) = fft_gavr.easy.stderr(idx,:)*zstar;
        civals(:,2) = fft_gavr.difficult.stderr(idx,:)*zstar;
        X=[freq(idx_freq);freq(idx_freq)]';

        clrs = [[.2 .2 1];[0 .7 0]];

        h = figure;
        set(h, 'Color', 'w');
        hh = plot_ci_filled(X, avrvals(idx_freq,:), avrvals(idx_freq,:)+civals(idx_freq,:), avrvals(idx_freq,:)-civals(idx_freq,:), clrs);

        set(hh,'LineWidth', 1.3);
        hh = title(sprintf('Spectral power over left lane change events [%s]', chan));
        set(hh,'FontSize', 24);
        hh = xlabel('Frequency (Hz)');
        set(hh,'FontSize', 20);
        hh = ylabel('Power');
        set(hh,'FontSize', 20);

        resize_window(h, [1000 600], [200 200]);

        hh = legend([{'Easy'},{'Difficult'}]);
        hh.FontSize = 20;

        ax = gca;
        ax.FontSize = 18;

    end

end