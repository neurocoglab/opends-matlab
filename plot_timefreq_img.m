function [ h ] = plot_timefreq_img( timelockstats, plot_params, todir, showfigs )

if nargin < 3 
    todir = [];
end

if nargin < 4
    showfigs = true;
end

time = timelockstats.time;
freq = timelockstats.freq;

channels = plot_params.channels;
transparency = plot_params.transparency;

dt = time(end)-time(1);
tdiv = dt / 20;
if tdiv+0.1>1, tdiv = ceil(tdiv); end
xtick = round(time(1)) : tdiv : round(time(end));
% xtick = plot_params.xtick;
ytick = plot_params.ytick;
ttitle = plot_params.title;

tticks = [];
for i = 1 : length(xtick)
    idx = find(time==xtick(i));
    if isempty(idx)
       [~,idx] = min(abs(time-xtick(i)));
    end
    tticks = [tticks idx];
end

fticks = [];
for i = 1 : length(ytick)
    idx = find(freq==ytick(i));
    if isempty(idx)
       [~,idx] = min(abs(freq-ytick(i)));
    end
    fticks = [fticks idx];
end

start_pos = [200 200];
idx_pos_clus = find([timelockstats.posclusters.prob]<0.05);
idx_neg_clus = find([timelockstats.negclusters.prob]<0.05);

for c = 1 : length(channels)

    idx_channel = find(strcmp(timelockstats.label,channels{c}),1);
    stat = squeeze(timelockstats.stat(idx_channel,:,:));
    
    mask = zeros(size(timelockstats.posclusterslabelmat,2),size(timelockstats.posclusterslabelmat,3));
    for j = 1 : length(idx_pos_clus)
        mask(squeeze(timelockstats.posclusterslabelmat(idx_channel,:,:)) == idx_pos_clus(j))=1;
    end
    for j = 1 : length(idx_neg_clus)
        mask(squeeze(timelockstats.negclusterslabelmat(idx_channel,:,:)) == idx_neg_clus(j))=1;
    end
    
    bmask = logical(mask>0);
    mask(mask==0) = transparency;
%     bmask = logical(mask>0); % squeeze(timelockstats.mask(idx_channel,:,:));
    
    if showfigs
        h = figure;
    else
        h = figure('Visible','off'); 
    end
    h.Color = 'w';
    hh = imagesc(stat);
    hh.AlphaData = mask;
    hold on;
    bw = bwperimtrace(bmask');
    for i = 1 : length(bw)
       plot(bw{i}(:,1),bw{i}(:,2),'Color',[130 38 31]/255);
       hold on;
    end
    
    ax = gca;
    
    ax.Title.String = sprintf('%s: [%s]', ttitle, channels{c});
    ax.Title.FontSize = 16;
    
    ax.YDir = 'normal';
    
    ax.XTick = tticks;
    ax.XTickLabel = xtick; %ttick_labels;
    ax.XLabel.String = 'Time (s)';
    ax.XLabel.FontSize = 14;
    
    ax.YTick = fticks;
    ax.YTickLabel = ytick; % ftick_labels;
    ax.YLabel.String = 'Freq (Hz)';
    ax.YLabel.FontSize = 14;
    
    caxis(plot_params.clim);
    if plot_params.colorbar
        colorbar;
    end
    
    resize_window(h, plot_params.window_size, start_pos);
    
    if max(start_pos) < 1000
        start_pos = start_pos + 10;
    end
    
    if ~isempty(todir)
        saveas(h,sprintf('%s/timefreq_%s_%s.png', todir, plot_params.name, channels{c}));
        if ~showfigs
            close(h);   
        end
    end

end

end

