%% Runs ICA on a subject and prompts for components to remove

addpath '/Users/lpzatr/OneDrive/OneDrive UoN/OneDrive - The University of Nottingham/synched/MATLAB/lib/fieldtrip-20180110'
addpath '/Users/lpzatr/OneDrive/OneDrive UoN/OneDrive - The University of Nottingham/synched/MATLAB/lib/cPCOH/functions';
ft_defaults

skip_existing = true;

% Get params
load processing_eeg_params_osx.mat
preproc = load('preproc_params_osx.mat');
proc = load('processing_params.mat');

load(proc.params.qc.file);
qc_score = cell2mat(qc_eeg(:,1));
idx = qc_score>1; %=proc.params.qc.cutoff;
subjects = qc_eeg(idx,2);

for s = 1 : length(subjects)
    subject = subjects{s};
    
    fprintf('\n\n== Starting subject %s ==\n\n', subject);
    outdir = sprintf('%s/%s/%s', preproc.params.root_dir, preproc.params.output_dir, subject);
    figdir = sprintf('%s/figures',outdir);
    
    flag_file = sprintf('%s/ica.done',outdir);
    results_file = sprintf('%s/processing_results_eeg_ica.mat',outdir);
    
    data = [];
    
    % Load data
    if exist(flag_file, 'file') && skip_existing
        fprintf('\nICA already done for %s.\n', subject);
    else
        data = load_eeg_data( params, preproc, subject );
    end
    
    % Perform ICA
    if isempty(data)
        fprintf('\nSkipping ICA for %s.\n', subject);
    else
        
        if exist(flag_file, 'file')
           delete(flag_file); 
        end
        
        cfg = params.eeg.cfg;
        cfg.dataset = sprintf('%s/%s/%s-eeg/%s.eeg', params.eeg.data_dir, subject, subject, subject);
        
        cfg.continuous = 'yes';
        cfg.hpfilter = 'yes';
        cfg.hpfreq = 1;
        cfg.hpfiltord = 6;
        cfg.lpfilter = 'yes';
        cfg.lpfreq = 60;
        cfg.lpfiltord = 6;
        cfg.demean = 'no';
        cfg.reref = 'yes';
        cfg.refchannel = 'all';

        data_ica = ft_preprocessing(cfg);
        
        cfg = [];
        cfg.method = 'runica'; % this is the default and uses the implementation from EEGLAB

        results.eeg.ica.comp = ft_componentanalysis(cfg, data_ica);

        cfg = [];
        cfg.colormap = 'jet';
        cfg.component = 1:20;       % specify the component(s) that should be plotted
        cfg.layout    = 'acticap-64ch-standard2.mat'; % specify the layout file that should be used for plotting
        cfg.comment   = 'no';
        h = figure;
        h.Color = 'w';
        ft_topoplotIC(cfg, results.eeg.ica.comp);
        resize_window(h,[1000 1000]);
        
        cfg.viewmode = 'component';
        ft_databrowser(cfg, results.eeg.ica.comp);
        resize_window(gcf,[1500 900]);
        
        % Prompt for bad components
        channels = input('Enter ICA components to remove: ','s');
        to_rem = cellfun(@str2num, strsplit(channels, ' '));

        if length(to_rem) > 0
           fprintf('Removing %d components...\n', length(to_rem));
           cfg = [];
           cfg.component = to_rem;
           data_flt = ft_rejectcomponent(cfg, results.eeg.ica.comp, data_ica);
           results.eeg.ica.removed = to_rem;
        end

        fprintf('Done. Saving results...\n');
        save(results_file, 'results', 'data_flt');
        saveas(gcf,sprintf('%s/ica_browser.fig',figdir));
        saveas(h,sprintf('%s/ica_topoplot.png',figdir));
        
        fid = fopen(flag_file,'w+');
        fclose(fid);

        close all;
        
    end
    
end

fprintf('\n\n===== DONE ALL SUBJECTS. HAVE A COFFEE ======\n');
