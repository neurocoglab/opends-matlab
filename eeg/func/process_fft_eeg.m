function results_fft = process_fft_eeg(cfg,data)

%% check config

% frequencies
if ~isfield(cfg,'freq')
    error('No frequencies specified')
end

% method
if ~isfield(cfg, 'method')
    cfg.method = 'fft';
end
if strcmp(cfg.method, 'fft') && length(cfg.freq)>2 %
    warning('Reducing freqs to foi')
    cfg.foi = [cfg.freq(1),cfg.freq(end)];
end

% 1/f correction?
if ~isfield(cfg, 'fcor')
    cfg.fcor = false;
end

% time selection?
if ~isfield(cfg,'toi') || isempty(cfg.toi)
    warning('No toi specified, computing fft for whole trial')
    cfg.toi = [];
elseif ~iscell(cfg.toi)
    cfg.toi = {cfg.toi};
end

% channel selection?
if ~isfield(cfg,'channels')
    cfg.channels = 1:length(data.label);
end

if ~isfield(cfg,'events')
    cfg.events = 1:size(data.event,1);
end
if length(cfg.toi) < length(cfg.events)
    error('toi is not specified for all events')
end

% trial selection?
if ~isfield(cfg,'trials')
    cfg.trials = 1:length(data.trial);
end

%% compute spectra

% prep output struct
results_fft = rmfield(data, 'trial');

% prep for loop
nevents = length(cfg.events);
ntrials = length(cfg.trials);
nchannels = length(cfg.channels);

for tr = 1:ntrials
    
    for ev = 1:nevents
        if isempty(data.event{ev,tr}) || isnan(data.event{ev,tr})
            continue
        end
        if ~isempty(cfg.toi{ev})
            tlog = data.time{cfg.trials(tr)}>=data.event{ev,tr}+cfg.toi{ev}(1) ...
                & data.time{cfg.trials(tr)}<=data.event{ev,tr}+cfg.toi{ev}(2);
        else
            tlog = ones(size(data.time{tr}));
        end
        
        for ch = 1:nchannels
            % ignore NaNs
            tid = tlog & ~isnan(data.trial{cfg.trials(tr)}(cfg.channels(ch),:));
            
            % compute spectrum
            switch cfg.method
                case 'fft'
                    dumfft = fft(data.trial{cfg.trials(tr)}(cfg.channels(ch),tid));
                    freqs = linspace(0,data.fsample,length(dumfft));
                    fid = find(freqs>=cfg.foi(1),1,'first'):find(freqs>=cfg.foi(2),1,'first');
                    datspctrm = abs(dumfft(fid));
                    freqs = freqs(fid);
                case 'pmtm'
                    %                 [dumfft,freqs] = pmtm(data.trial{cfg.trials(tr)}(cfg.channels(ch),tid),2,...
                    %                     sum(tid),data.fsample);
                    %                 fid = find(freqs>=cfg.foi(1),1,'first'):find(freqs>=cfg.foi(2),1,'first');
                    %                 datspctrm = dumfft(fid);
                    %                 freqs = freqs(fid);
                    freqs = cfg.freq;
                    datspctrm = pmtm(data.trial{cfg.trials(tr)}(cfg.channels(ch),tid),5,...
                        freqs,data.fsample);
            end
            
            %         figure;
            %         plot(freqs,datspctrm)
            
            % do 1/f correction
            if cfg.fcor
                % fit 1/f - is this robust?
                p = polyfit(log10(freqs), reshape(log10(datspctrm),size(freqs)), 1);
                
                % remove 1/f
                datspctrm = datspctrm - reshape(10.^(log10(freqs)*p(1)+p(2)),size(datspctrm));
                
                %             hold on; plot(freqs,reshape(10.^(log10(freqs)*p(1)+p(2)),size(datspctrm)),'linewidth',2)
            end
            
            % store data
            results_fft.powspctrum(ev,tr,ch,:) = datspctrm;
%             results_fft.powspctrm{tr}(ch,:) = datspctrm;
        end
    end
end

%% admin
results_fft.freq = freqs;
cfg.freq = freqs;
results_fft.cfg_fft = cfg;

end