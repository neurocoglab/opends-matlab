function [ h ] = plot_luminance( results, data, params, out2file )

if nargin < 4
    out2file = false;
end

 % Output plots to images
if out2file
    h = figure('visible','off');
    outdir = sprintf('%s/%s/%s/%s/figures', params.root_dir, params.data_dir, params.output_dir, data.subject);   
    if ~exist(outdir,'dir')
       mkdir(outdir); 
    end
else
    h = figure;
end
h.Color = 'w';

idx_offset = params.luminance.use_offset;
if idx_offset < 1
   [~,idx_offset] = max(results.luminance.r2); 
end

t_offset = params.luminance.offsets(idx_offset);
T_s = 1 / data.eye.Fs;
N_offset = length(params.luminance.offsets);

% Remove outliers and smooth if necessary
idx_keep = find(abs(zscore(results.diam_left(results.luminance.idx))) < params.luminance.outlier_lim);
X = results.luminance.y(idx_keep);
y = results.diam_left(idx_keep);

offset = round(t_offset / T_s);
% Shift luminance by offset
X = X(max(offset+1,1):min(end+offset,length(X)));
y = y(1:length(X));

%if params.luminance.downsample > 1
    % Downsample by factor (prevents rank deficiency)
    idx_ds = 1:round(length(X)/2000):length(X);
    X = X(idx_ds);
    y = y(idx_ds);
    
%end

y = y(~isnan(X));
X = X(~isnan(X));
if params.luminance.smooth > 0
    X = smooth(X, params.luminance.smooth);
    y = smooth(y, params.luminance.smooth);
end

scatter(zscore(X),zscore(y));
hh = lsline;
hh.Color = 'r';
% xlim([-4 4]);
ylim([-4 4]);
box on;
hh = xlabel('Luminance (z-score)');
hh.FontSize = 12;
hh = ylabel('Pupil diameter (z-score)');
hh.FontSize = 12;
hh = title(sprintf('Luminance effect (%1.1fs offset) - Subject %s', t_offset,data.subject));
hh.FontSize = 16;
if out2file
    saveas(h, sprintf('%s/lum_scatter.fig', outdir));
    saveas(h, sprintf('%s/lum_scatter.png', outdir));
    close(h);
end

if out2file
    h = figure('visible','off');
else
    h = figure;
end
h.Color = 'w';
tmin = results.luminance.ts{idx_offset}/1000/60;
plot(tmin,zscore(results.diam_left(results.luminance.idx_offset{idx_offset})));
hold on; plot(tmin,zscore(results.luminance.diam_left{idx_offset}));
hh = xlabel('Time (min)');
hh.FontSize = 12;
hh = ylabel('Pupil diameter (z-score)');
hh.FontSize = 12;
hh = title(sprintf('Luminance correction (%1.1f s offset) - Subject %s', t_offset,data.subject));
hh.FontSize = 16;
legend([{'Original'},{'Corrected'}]);
ylim([-3 4.5]);
% xlim([0 3]);
resize_window(h,[1000,500]);
if out2file
    saveas(h, sprintf('%s/lum_corrected_pd.fig', outdir));
    xlim([0 3]);
    saveas(h, sprintf('%s/lum_corrected_pd.png', outdir));
    close(h);
end


% Plot R2 over offsets (cross-regression)

N_offset = length(results.luminance.lm);
R2 = zeros(N_offset,1);

for i = 1 : N_offset
    R2(i) = results.luminance.lm{i}.Rsquared.Ordinary;
end

% R2 = results.luminance.r2;

if out2file
    h = figure('visible','off');
else
    h = figure;
end
h.Color = 'w';

yyaxis left;
plot(params.luminance.offsets, results.luminance.r2);
hh = ylabel('R^2');
hh.FontSize = 12;
hold on;
yyaxis right;
plot(params.luminance.offsets, results.luminance.rmse);
hh = ylabel('RMSE');
hh.FontSize = 12;

% R^2 versus Offsets 
hh = xlabel('Offset (s)');
hh.FontSize = 12;

hh = title(sprintf('Regression (Luminance v. PD) - Subject %s', data.subject));
hh.FontSize = 16;
% ylim([0 0.3]);
resize_window(h,[1000,500]);
if out2file
    saveas(h, sprintf('%s/lum_r2_offsets.fig', outdir));
    saveas(h, sprintf('%s/lum_r2_offsets.png', outdir));
    close(h);
end


end

