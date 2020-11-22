%{ 
select_data_eeg 

Inputs:
- cfg
- dataIn: FieldTrip structure of either field data or power spectra

The cfg specifies characteristics of desired trials and/or time windows
around events. 
To select time windows, specify:
- cfg.event: (required) The event to keep (only one event can be used)
- cfg.toi: (required) Time points around this event to keep, either as:
  a time series (array) or as [time_min, time_max].
  Note that times are in seconds relative to the event!
- cfg.dt: if specifying toi as [time_min, time_max], you can also specify
  a separate time step dt. If dt is not given, the original dt of the data
  is used

Output: 
- dataOut: FieldTrip struct with only desired data

%}

function [dataOut,keepID] = select_data_eeg(cfg,dataIn)

% check settings
if isfield(cfg, 'toi') && ~isfield(cfg, 'event') || ...
        ~isfield(cfg, 'toi') && isfield(cfg, 'event')
    error('Event and toi need to both be specified to select data')
end

% find datatype
if isfield(dataIn, 'powspctrm')
    datatype = 'spectrum';
else
    datatype = 'field';
end


%% select rounds

keepID = true(1,size(dataIn.time,2)); % start with all trials

% only use trials specified in config
if isfield(cfg,'trial')
    keepID(setdiff(find(keepID),cfg.trials)) = false;
end

% only use specified rounds
if isfield(cfg,'rounds')
    keepID(dataIn.rounds~=cfg.rounds) = false;
end

% only use specified repeats
if isfield(cfg,'repeats')
    keepID(dataIn.repeats~=cfg.repeats) = false;
end

% trial admin

% create output struct
dataOut = dataIn;
fields = fieldnames(dataOut);
for f = 1:length(fields)
    if eval(['size(dataIn.',fields{f}, ',2) == length(keepID)'])
        eval(['dataOut.',fields{f}, ' = dataIn.',fields{f},'(:,keepID);'])
    end
end

%% select time wimdows around events within trial

if isfield(cfg, 'toi') && isfield(cfg, 'event')
    
    rmfields = {'time','event','event_info'};
    dataOutTmp = rmfield(dataOut,rmfields);
    
    fields = fieldnames(dataOut);
    fields = fields(~ismember(fields,rmfields));
    
    % find time and index step
    if length(cfg.toi) == 2 && ~isfield(cfg, 'dt')
        % if dt is not specified, keep the original time step
        cfg.dt = 1/dataOut.fsample;
        cfg.toi = cfg.toi(1):cfg.dt:cfg.toi(2);
    elseif length(cfg.toi) > 2
        % update dt
        cfg.dt = cfg.toi(2)-cfg.toi(1);
    end
    DT = round(cfg.dt*dataOut.fsample); % index stepsize
    
    % select relevant data for each trial
    ctr = 1;
    for tr = 1:length(dataOut.time)
        
        % in case there are several events of this type in the trial
        for ev = 1:length(dataOut.event{cfg.event,tr})
        
            for f = 1:length(fields)
                if eval(['size(dataOutTmp.',fields{f}, ',2) == size(dataOut.time,2)'])
                    eval(['dataOutTmp.',fields{f}, '(:,ctr) = dataOut.',fields{f},'(:,tr);'])
                end
            end
            
            % find event ID
            [~,eventID] = min(abs(dataOut.time{tr}-dataOut.event{cfg.event,tr}(ev)));
        
            % find indices of time window around event
            selectIDs = eventID+round(cfg.toi(1)*dataOut.fsample): DT: eventID+round(cfg.toi(end)*dataOut.fsample);
            
            if sum(selectIDs < 1)
                error('Window starts too early')
            elseif selectIDs(end) > length(dataOut.time{tr})
                error('Window ends too late')
            end
            
            if isfield(dataIn, 'trial')
                dataOutTmp.trial{ctr} = dataOut.trial{tr}(:,selectIDs);
            elseif isfield(dataIn, 'powspctrm')
                dataOutTmp.powspctrm{ctr} = dataOut.powspctrm{tr}(:,:,selectIDs);
            end
            
            % admin
            % set time axis
            dataOutTmp.time{ctr} = dataOut.time{tr}(selectIDs) - dataOut.time{tr}(eventID);
            dataOutTmp.event{cfg.event,ctr} = 0;
            
            % TODO: this should be changed!
            dataOutTmp.event_info{1,ctr} = dataOut.event_info{cfg.event,tr}(ev);
            
        
        ctr = ctr +1;
        end
    end
   
    clear dataOut
    dataOut = dataOutTmp;

end

end
