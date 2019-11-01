function [ h ] = plot_eeg_erp_stats ( stats, params, figdir, show_plots )
%plot_eeg_erp_stats 

%     '#1f77b4',  // muted blue
%     '#ff7f0e',  // safety orange
%     '#2ca02c',  // cooked asparagus green
%     '#d62728',  // brick red
%     '#9467bd',  // muted purple
%     '#8c564b',  // chestnut brown
%     '#e377c2',  // raspberry yogurt pink
%     '#7f7f7f',  // middle gray
%     '#bcbd22',  // curry yellow-green
%     '#17becf'   // blue-teal

plotly_clrs = [{'#1f77b4'},{'#ff7f0e'},{'#2ca02c'}];

fields = fieldnames(stats);

if nargin < 3
   figdir = []; 
end

if nargin < 4
    show_plots = true;
end

for f = 1 : length(fields)
    
    ff = fields{f};
    params_f = params.eeg.erp.plots.(ff);
    params_f.varname = ff;
    
    plot_topo(stats.(ff), params_f, ff);
    
    params_f.zlim = params_f.zlim_comp;
    
%     if ~isempty(figdir), plot_erp(stats.(ff), params_f, {ff}); end
    
    % Plot any comparisons
    subfields = fieldnames(stats.(ff));
    subfields(strcmp(subfields,'timelockstats'))=[];
    subfields(strcmp(subfields,'timelock_avr'))=[];
    
    for l = 1 : length(subfields)
        ll = subfields{l};
        levels = stats.(ff).(ll).levels;
        plot_topo(stats.(ff).(ll), params_f, [ff '_' ll]);
        params_f.varname = [ff '_' ll];
%         if ~isempty(figdir), plot_erp(stats.(ff).(ll), params_f, stats.(ff).(ll).levels); end
    end
    
end

    % Plots the topology of the t-statistic over time windows
    function h = plot_topo ( statf, params_f, varname )

        stat = statf.timelockstats;
        tlock_ga = statf.timelock_avr;

        if length(tlock_ga) == 2
            cfg_m = [];
            cfg_m.operation = 'subtract';
            cfg_m.parameter = 'avg';
            tlock_ga = ft_math(cfg_m, tlock_ga{1}, tlock_ga{2});
        end

        alpha = 0.025;
        pos = false(size(stat.prob)); neg = false(size(stat.prob));
        
        if ~isempty(stat.posclusters)
            pos_cluster_pvals = [stat.posclusters(:).prob]; 
            pos_signif_clust = find(pos_cluster_pvals < alpha);
            pos = ismember(stat.posclusterslabelmat, pos_signif_clust);
        end
        if ~isempty(stat.negclusters)
            neg_cluster_pvals = [stat.negclusters(:).prob];
            neg_signif_clust = find(neg_cluster_pvals < alpha);
            neg = ismember(stat.negclusterslabelmat, neg_signif_clust);
        end

        tmin = params_f.tlim(1); tmax = params_f.tlim(2);
        timestep = params_f.tstep; % timestep between time windows for each subplot (in seconds)
        sampling_rate = params.eeg.erp.resamplefs; 
        j = tmin:timestep:tmax;
        m = round(find(stat.time>tmin,1):timestep*sampling_rate:find(stat.time>tmax,1)-1);
        if length(m) < length(j) 
           m = round([m m(end)+timestep*sampling_rate]);
        end

        [i1,i2] = match_str(tlock_ga.label, stat.label);

        Nj = length(j)-1;
        Ncol = round(sqrt(Nj));
        Nrow = ceil(Nj/Ncol);

        if ~show_plots
            h = figure('Visible', 'off');
        else
            h = figure;
        end
        h.Color = 'w';


        for k = 1 : Nj
           subplot(Nrow,Ncol,k);
           cfg = [];
           cfg.xlim=[j(k) j(k+1)];   % time interval of the subplot
           cfg.zlim = params_f.zlim;
           % If a channel reaches this significance, then
           % the element of pos_int with an index equal to that channel
           % number will be set to 1 (otherwise 0).

           % Next, check which channels are significant over the
           % entire time interval of interest.
           pos_int = zeros(numel(tlock_ga.label),1);
           neg_int = zeros(numel(tlock_ga.label),1);
           pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
           neg_int(i1) = all(neg(i2, m(k):m(k+1)), 2);

           cfg.highlight = 'on';
           % Get the index of each significant channel
           cfg.highlightchannel = find(pos_int | neg_int);
           cfg.highlightsymbol = '.';
           cfg.highlightsize = 14;

           cfg.comment = 'xlim';
           cfg.commentpos = 'title';
%            cfg.layout = 'acticap-64ch-standard2.mat';
           cfg.layout = 'biosemi64.lay';
           cfg.interactive = 'no';
           cfg.colorbar = 'no';
           cfg.parameter='stat';
           
           evalc('ft_topoplotER(cfg, stat);');

%            evalc('ft_topoplotER(cfg, tlock_ga);');

           colormap(brewermap(256, '*RdYlBu'));

           hh = gca;
           hh.Title.String = sprintf('%d to %d ms', round(j(k)*1000), round(j(k+1)*1000));
           hh.Title.FontSize = 15;
        end

        resize_window(h, [1000 1000], [400 400]);

        hh = colorbar;
        hh.Position = [.93 0.4 .015 .3];
        hh.Ticks = [cfg.zlim(1) cfg.zlim(1)+diff(cfg.zlim)/2,cfg.zlim(2)];
        hh.TickLabels = compose('%1.1f',hh.Ticks');  %[{'-4.0'},{'0.0'},{'4.0'}];
        hh.FontSize = 18;

        if ~isempty(figdir)
            fig_file = sprintf('%s/erp_topo_%s.png', figdir, varname);
            saveas(h, fig_file);
            fprintf(' Saved %s topo plot to %s\n', varname, fig_file);
        end
        
        tt = params_f.title;
        
        if isfield(statf, 'levels')
            tt = sprintf('%s (%s - %s)', tt, statf.levels{1}, statf.levels{2});
        end
        
        hh = suptitle(tt);
        hh.FontSize = 24;
        hh.FontWeight = 'bold';

        if ~show_plots
            close(h);
        end


    end

    % Plots the ERP as lines with confidence intervals
    function h = plot_erp ( statf, params_f, groups )
        
        t_win = params_f.tlim;
%         td = diff(t_win);
%         t_win(1) = t_win(1) - td*1.2;
%         t_win(2) = t_win(2) + td*1.2;
        
        if iscell(statf.timelock_avr)
            N_subj = statf.timelock_avr{1}.dof(1,1);
            ttime = statf.timelock_avr{1}.time;
        else
            N_subj = statf.timelock_avr.dof(1,1);
            ttime = statf.timelock_avr.time;
        end
        
        idx1 = find(ttime < t_win(1), 1, 'last');
        idx2 = find(ttime > t_win(2), 1);
        if isempty(idx1), idx1 = 1; end
        if isempty(idx2), idx2 = length(ttime); end
        
        for c = 1 : length(params_f.erp_channels)
            
            ch = params_f.erp_channels{c};
            idx_c = find(strcmp(statf.timelockstats.label, ch));
            
            s = [];
            s.mask = statf.timelockstats.mask(idx_c,:);
            idxpos = find(([statf.timelockstats.posclusters.prob] > 0.025) & ([statf.timelockstats.posclusters.prob] < 0.05));
            idxneg = find(([statf.timelockstats.negclusters.prob] > 0.025) & ([statf.timelockstats.negclusters.prob] < 0.05)); 
            s.mask2 = ismember(statf.timelockstats.posclusterslabelmat(idx_c,:), idxpos) | ...
                      ismember(statf.timelockstats.negclusterslabelmat(idx_c,:), idxneg);
            s.mask2(s.mask) = false;
            s.time = statf.timelockstats.time;
            stats_c = [];
            stats_c.timelockstats = s;
            stats_c.timelock_avr = [];
            if iscell(statf.timelock_avr)
                for i = 1 : length(statf.timelock_avr)
                    s = [];
                    s.avg = statf.timelock_avr{i}.avg(idx_c,idx1:idx2);
                    s.var = statf.timelock_avr{i}.var(idx_c,idx1:idx2);
                    s.time = statf.timelock_avr{i}.time(idx1:idx2);
                    stats_c.timelock_avr = [stats_c.timelock_avr {s}];
                end
            else
                s = [];
                s.avg = statf.timelock_avr.avg(idx_c,idx1:idx2);
                s.var = statf.timelock_avr.var(idx_c,idx1:idx2);
                s.time = statf.timelock_avr.time(idx1:idx2);
                stats_c.timelock_avr = {s};
            end
            
            sigdims = [-4 1];
%             sigdims = [-10 50];
            sig_clrs = [{['#2ca02c' dec2hex(40)]},{['#d62728' dec2hex(20)]}];
            h = plot_timelocked(stats_c, N_subj, groups, params.eeg.erp.plots.layouts.boxplot, sigdims, sig_clrs);
            h.layout.title = sprintf('%s - %s', params_f.title, ch);
            h.layout.yaxis.range = [-5 4]; % [min(y(:)) max(y(:))]*1.3;
            if length(groups) > 1
                h.layout.showlegend = true;
                h.layout.legend = struct('x', 0.8, 'y', 0.1); 
            end
            
            h.PlotOptions.FileName = sprintf('%s/erp_%s_%s', figdir, params_f.varname, ch);
            plotlyoffline(h);

            if show_plots
                web(['file://' h.PlotOptions.FileName '.html'], '-new', '-notoolbar');
            end

            % Gives a URL error $%#^!
            % Workaround is to save manually...
%             saveplotlyfig(h, sprintf('%s/erp_%s.png', figdir));
          
        end
        
        
    end



end