%% Load EEG markers and match to Simulation Log events

clear;

addpath('../lib/fieldtrip-20190419/');
ft_defaults;

load preproc_params_osx.mat;

% Temp DEBUG
params.clobber=true;

% subjects = csvimport(params.subject_list_file, 'noHeader', true, 'outputAsChar', true);

data_dir = '/Users/lpzatr/Documents/data/driving';
sim_dir = sprintf('%s/processed', data_dir);
eeg_dir = sprintf('%s/eeg', data_dir);
subjects = [{'0739'},{'1234'}];

for i = 1 : length(subjects)
    
    subject = subjects{i};
    
    % Load in marker file
    cfg = [];
    cfg.dataset = sprintf('%s/%s/%s-eeg/%s.eeg', eeg_dir, subject, subject, subject);
	eeg_data = ft_preprocessing(cfg);
    mrk_file = sprintf('%s/%s/%s-eeg/%s.vmrk', eeg_dir, subject, subject, subject);
    eeg_events = ft_read_event(mrk_file);
    idx_stim = strcmp({eeg_events.type},'Stimulus');
    eeg_events = eeg_events(idx_stim);
    t_eeg = eeg_data.time{1}([eeg_events.sample]);
    ft_write_event
    
    % Load in simulation data
    sim_file = sprintf('%s/%s/%s_events-All.csv', sim_dir, subject, subject);
     
    exist(sim_file, 'file')
    events_table = readtable(sim_file);
    t_sim = events_table.Millis;
    t_start = events_table.Millis(find(strcmp(events_table.EventType,'SimulatorStarted')),:);
    t_sim = t_sim - t_start;
    idx_lane = find(strcmp(events_table.EventType,'LaneChange'));
    t_lane = t_sim(idx_lane);
    
    events_file = sprintf('%s/%s/%s_processing_results.mat', sim_dir, subject, subject);
    load(events_file);
    
    for j = 1 : length(results.events.left_change.events)
        
        event_j = results.events.left_change.events(j);
        idx = find(abs(t_lane - event_j) < 1000);
        fprintf('%1.4f\n', t_lane(idx(1)) - event_j);
        
    end
    
    
end