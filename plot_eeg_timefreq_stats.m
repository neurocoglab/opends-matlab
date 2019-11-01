function [ h ] = plot_eeg_timefreq_stats ( stats, params, figdir, show_plots )

if nargin < 3
   figdir = []; 
end

if nargin < 4
    show_plots = true;
end

fields = fieldnames(stats);

for f = 1 : length(fields)
    
    var_f = fields{f};
    
    params_f = params.eeg.timefreq.plots.(var_f);
    params_f.name = [params_f.name '_all'];
    h = plot_timefreq_img(stats.(var_f).timelockstats, params_f, figdir, show_plots);
    params_f.name = [params_f.name '_diff'];
    h = plot_timefreq_img(stats.(var_f).diff.timelockstats, params_f, figdir, show_plots);
    params_f.name = [params_f.name '_outcome'];
    h = plot_timefreq_img(stats.(var_f).outcome.timelockstats, params_f, figdir, show_plots);
    
    
end





end