function [ h ] = plot_blinks( results, data, params, out2file )

if nargin < 4
    out2file = false;
end

time_sec = results.t / 60000;
            
gap_clr = [0.9 0.8 0.8];
blink_clr = [0.8 0.9 0.8];

if out2file
    h = figure('visible','off');
else
    h = figure;
end
set(h, 'Color', 'w');

% Blinks
for j = 1 : length(results.blink_ints)
   int_j = results.blink_ints{j};

   for k = 1 : size(int_j,1)
       x1 = time_sec(int_j(k,1));
       try
       x2 = time_sec(int_j(k,1)+int_j(k,2));
       catch
           a=0;
       end
       hr = rectangle('Position',[x1 -500 x2-x1 1000], ...
                      'EdgeColor','w', ...
                      'FaceColor', blink_clr); 
   end
end

% Gaps
for j = 1 : size(data.eye.tgap,1)
   x1 = time_sec(data.eye.tgap(j,1));
   x2 = x1 + data.eye.tgap(j,2) / 60000;
   hr = rectangle('Position',[x1 -500 x2-x1 1000], ...
                  'EdgeColor','w', ...
                  'FaceColor', gap_clr); 
end

hold on;
hh = plot(time_sec, data.eye.diam_left);
set(hh,'Color',[0.7 0.7 0.7]);

plot(time_sec, results.diam_left, 'b');

% Derivative + thresholds
offset = 10; if isfield(params.blink,'plot_offset'), offset=params.blink.plot_offset; end
plot(time_sec, results.d_x*2 + offset, 'm');

thr_color = [0.5 0.5 0.5];
thres1 = params.blink.thres(1)*2 + offset;
line([time_sec(1) time_sec(end)], [thres1 thres1], ...
     'Color',thr_color, 'LineStyle', '--');
thres1 = params.blink.thres(2)*2 + offset;
line([time_sec(1) time_sec(end)], [thres1 thres1], ...
     'Color',thr_color, 'LineStyle', '--');

hh = title('Pupil diameter blink corrected');
set (hh, 'FontSize', 14);
ylim([min(results.diam_left),max(results.diam_left)]);
h.Position = [400 400 1500 500];
% resize_window(h, [1500 500], [400 400]);

if out2file
    outdir = sprintf('%s/%s/%s/%s/figures', params.root_dir, params.data_dir, params.output_dir, data.subject);
    if ~exist(outdir,'dir')
       mkdir(outdir); 
    end
    saveas(h, sprintf('%s/blinks.fig', outdir));
    xlim([0 3]);
    saveas(h, sprintf('%s/blinks.png', outdir));
    close(h);
end

end

