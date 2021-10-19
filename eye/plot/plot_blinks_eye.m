function [ h ] = plot_blinks_eye( data, params, out2file )

if nargin < 3
    out2file = false;
end

time_sec = data.eye.t / 60000;
            
gap_clr = params.eye.gaps.plots.color;
blink_clr = params.eye.blinks.plots.color;

if out2file
    h = figure('visible','off');
else
    h = figure;
end
set(h, 'Color', 'w');

% Blinks
for j = 1 : size(data.eye.blinks.intervals,1)
   x1 = time_sec(data.eye.blinks.intervals(j,1));
   try
   x2 = time_sec(data.eye.blinks.intervals(j,2));
   catch
       a=0;
   end
   hr = rectangle('Position',[x1 -500 x2-x1 1000], ...
                  'EdgeColor','w', ...
                  'FaceColor', blink_clr); 
end

% for j = 1 : length(data.eye.blinks.blink_ints)
%    int_j = data.eye.blinks.blink_ints{j};
% 
%    for k = 1 : size(int_j,1)
%        x1 = time_sec(int_j(k,1));
%        try
%        x2 = time_sec(int_j(k,1)+int_j(k,2));
%        catch
%            a=0;
%        end
%        hr = rectangle('Position',[x1 -500 x2-x1 1000], ...
%                       'EdgeColor','w', ...
%                       'FaceColor', blink_clr); 
%    end
% end

% Gaps
for j = 1 : size(data.eye.tgap,1)
   x1 = time_sec(data.eye.tgap(j,1));
   x2 = x1 + data.eye.tgap(j,2) / 60000;
   hr = rectangle('Position',[x1 -500 x2-x1 1000], ...
                  'EdgeColor','w', ...
                  'FaceColor', gap_clr); 
end

hold on;
hh = plot(time_sec, data.eye.diam);
set(hh,'Color',[0.7 0.7 0.7]);

plot(time_sec, data.eye.blinks.diam, 'b');

% Derivative + thresholds
offset = 10; if isfield(params.eye.blinks.plots,'offset'), offset=params.eye.blinks.plots.offset; end
plot(time_sec, data.eye.blinks.d_x*2 + offset, 'm');

thr_color = [0.5 0.5 0.5];
thres1 = params.eye.blinks.thres(1)*2 + offset;
line([time_sec(1) time_sec(end)], [thres1 thres1], ...
     'Color',thr_color, 'LineStyle', '--');
thres1 = params.eye.blinks.thres(2)*2 + offset;
line([time_sec(1) time_sec(end)], [thres1 thres1], ...
     'Color',thr_color, 'LineStyle', '--');

hh = title('Pupil diameter blink corrected');
set (hh, 'FontSize', 14);
ylim([min(data.eye.diam),max(data.eye.diam)]);
h.Position = [400 400 1500 500];

xlabel('Time (min)');

if out2file
    outdir = sprintf('%s/%s/figures', params.io.output_dir, data.subject);
    if ~exist(outdir,'dir')
       mkdir(outdir); 
    end
    saveas(h, sprintf('%s/eye_blinks.fig', outdir));
    xlim([0 3]);
    saveas(h, sprintf('%s/eye_blinks.png', outdir));
    close(h);
end

end

