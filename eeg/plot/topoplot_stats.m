function [ h ] = topoplot_stats( cfg, channels, stat, stat_mask )
%TOPOPLOT_STATS Plots a given statistic as a layout topoplot

if nargin < 4
    stat_mask = [];
end

ft_defaults;

if ft_platform_supports('griddata-v4')
  default_interpmethod = 'v4';
else
  % Octave does not support 'v4', and 'cubic' is not yet implemented
  default_interpmethod = 'linear';
end

% Defaults
cfg.xlim              = ft_getopt(cfg, 'xlim',             'maxmin');
cfg.ylim              = ft_getopt(cfg, 'ylim',             'maxmin');
cfg.zlim              = ft_getopt(cfg, 'zlim',             'maxmin');
cfg.style             = ft_getopt(cfg, 'style',            'both');
cfg.gridscale         = ft_getopt(cfg, 'gridscale',         67);
cfg.interplimits      = ft_getopt(cfg, 'interplimits',     'head');
cfg.interpolation     = ft_getopt(cfg, 'interpolation',     default_interpmethod);
cfg.box               = ft_getopt(cfg, 'box',               'off');
cfg.contournum        = ft_getopt(cfg, 'contournum',        6);
cfg.colormap          = ft_getopt(cfg, 'colormap',         'default');
cfg.colorbar          = ft_getopt(cfg, 'colorbar',         'no');
cfg.colorbartext      = ft_getopt(cfg, 'colorbartext',    '');
cfg.shading           = ft_getopt(cfg, 'shading',          'flat');
cfg.comment           = ft_getopt(cfg, 'comment',          'auto');
cfg.fontsize          = ft_getopt(cfg, 'fontsize',          8);
cfg.fontweight        = ft_getopt(cfg, 'fontweight',       'normal');
cfg.baseline          = ft_getopt(cfg, 'baseline',         'no'); % to avoid warning in timelock/freqbaseline
cfg.trials            = ft_getopt(cfg, 'trials',           'all', 1);
cfg.interactive       = ft_getopt(cfg, 'interactive',      'no');
cfg.hotkeys           = ft_getopt(cfg, 'hotkeys',          'no');
cfg.renderer          = ft_getopt(cfg, 'renderer',          []); % let MATLAB decide on the default
cfg.marker            = ft_getopt(cfg, 'marker',           'on');
cfg.markersymbol      = ft_getopt(cfg, 'markersymbol',     'o');
cfg.markercolor       = ft_getopt(cfg, 'markercolor',       [0 0 0]);
cfg.markersize        = ft_getopt(cfg, 'markersize',        2);
cfg.markerfontsize    = ft_getopt(cfg, 'markerfontsize',    8);
cfg.highlight         = ft_getopt(cfg, 'highlight',        'off');
cfg.highlightchannel  = ft_getopt(cfg, 'highlightchannel', 'all', 1); % highlight may be 'on', making highlightchannel {} meaningful
cfg.highlightsymbol   = ft_getopt(cfg, 'highlightsymbol',  '*');
cfg.highlightcolor    = ft_getopt(cfg, 'highlightcolor',    [0 0 0]);
cfg.highlightsize     = ft_getopt(cfg, 'highlightsize',     6);
cfg.highlightfontsize = ft_getopt(cfg, 'highlightfontsize', 8);
cfg.labeloffset       = ft_getopt(cfg, 'labeloffset',       0.005);
cfg.maskparameter     = ft_getopt(cfg, 'maskparameter',     []);
cfg.component         = ft_getopt(cfg, 'component',         []);
cfg.directionality    = ft_getopt(cfg, 'directionality',    []);
cfg.channel           = ft_getopt(cfg, 'channel',          'all');
cfg.refchannel        = ft_getopt(cfg, 'refchannel',        []);
cfg.figurename        = ft_getopt(cfg, 'figurename',        []);
cfg.interpolatenan    = ft_getopt(cfg, 'interpolatenan',   'yes');
cfg.commentpos        = ft_getopt(cfg, 'commentpos',       'layout');
cfg.scalepos          = ft_getopt(cfg, 'scalepos',         'layout');
cfg.figure            = ft_getopt(cfg, 'figure',           'yes');

cfg.label = channels;
labels = cfg.layout.label;

stat_to_plot = nan(length(labels),1);
for i = 1 : length(labels)
    idx = find(strcmp(channels, labels{i}));
    if ~isempty(idx)
        stat_to_plot(i) = stat(idx);
    end
end

if strcmp(cfg.interplimits, 'head')
    interplimits = 'mask';
else
    interplimits = cfg.interplimits;
end

% Get physical min/max range of z:
if strcmp(cfg.zlim, 'maxmin')
    zmin = min(stat);
    zmax = max(stat);
elseif strcmp(cfg.zlim, 'maxabs')
    zmin = -max(max(abs(stat)));
    zmax = max(max(abs(stat)));
elseif strcmp(cfg.zlim, 'zeromax')
    zmin = 0;
    zmax = max(stat);
elseif strcmp(cfg.zlim, 'minzero')
    zmin = min(stat);
    zmax = 0;
else
    zmin = cfg.zlim(1);
    zmax = cfg.zlim(2);
end
  
opt = {};
opt = ft_setopt(opt, 'interpmethod',  cfg.interpolation);
opt = ft_setopt(opt, 'interplim',     interplimits);
opt = ft_setopt(opt, 'gridscale',     cfg.gridscale);
opt = ft_setopt(opt, 'outline',       cfg.layout.outline);
opt = ft_setopt(opt, 'shading',       cfg.shading);
opt = ft_setopt(opt, 'isolines',      cfg.contournum);
opt = ft_setopt(opt, 'mask',          cfg.layout.mask);
opt = ft_setopt(opt, 'style',         'surfiso');
opt = ft_setopt(opt, 'datmask',       stat_mask);
opt = ft_setopt(opt, 'box',           cfg.box);
opt = ft_setopt(opt, 'clim',          [zmin zmax]);

[Zi, h] = ft_plot_topo(cfg.layout.pos(:,1), cfg.layout.pos(:,2), stat_to_plot, opt{:});
ft_colormap(cfg.colormap);

switch cfg.marker
    case {'off', 'no'}
      % do not show the markers
    case {'on', 'labels', 'numbers'}
      channelsToMark = 1:length(stat);
      channelsToHighlight = [];
      for icell = 1:length(cfg.highlight)
        if ~strcmp(cfg.highlight, 'off')
          channelsToHighlight = [channelsToHighlight; match_str(channels, cfg.highlightchannel{icell})];
        end
      end
      if strcmp(cfg.interpolatenan, 'no')
        channelsNotMark = channelsToHighlight;
      else
        channelsNotMark = union(find(isnan(stat)), channelsToHighlight);
      end
      channelsToMark(channelsNotMark) = [];
      [dum, layoutindex] = match_str(ft_channelselection(channelsToMark, channels), cfg.layout.label);
      templay = [];
      templay.outline = cfg.layout.outline;
      templay.mask    = cfg.layout.mask;
      templay.pos     = cfg.layout.pos(layoutindex,:);
      templay.width   = cfg.layout.width(layoutindex);
      templay.height  = cfg.layout.height(layoutindex);
      templay.label   = cfg.layout.label(layoutindex);
      if strcmp(cfg.marker, 'labels') || strcmp(cfg.marker, 'numbers')
        labelflg = 1;
      else
        labelflg = 0;
      end
      if strcmp(cfg.marker, 'numbers')
        for ichan = 1:length(layoutindex)
          templay.label{ichan} = num2str(match_str(channels,templay.label{ichan}));
        end
      end
      ft_plot_layout(templay, 'box', 'no', 'label', labelflg, ...
        'pointsymbol',  cfg.markersymbol, ...
        'pointcolor',   cfg.markercolor, ...
        'pointsize',    cfg.markersize, ...
        'fontsize',     cfg.markerfontsize, ...
        'labeloffset',  cfg.labeloffset, ...
        'labelalignh', 'left', ...
        'labelalignv', 'bottom');
    otherwise
      ft_error('incorrect value for cfg.marker');
end

axis off;

% Plot colorbar
if isfield(cfg, 'colorbar')
if strcmp(cfg.colorbar, 'yes')
  c = colorbar;
  ylabel(c, cfg.colorbartext);
elseif ~strcmp(cfg.colorbar, 'no')
  c = colorbar('location', cfg.colorbar);
  ylabel(c, cfg.colorbartext);
end
end