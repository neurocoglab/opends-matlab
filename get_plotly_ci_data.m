function [ xy_line, xy_ci ] = get_plotly_ci_data( x, y, upper, lower, colours, name, alpha )
%get_plotly_ci Returns two plotly dataset for rendering as a shaded
%              confidence interval, defined by ci

if nargin < 7
   alpha = 125;
else
   alpha = round(255*alpha);
end

if nargin < 6
    colours = [{'r'},{'g'},{'b'}];
end

if nargin < 5
   name = ''; 
end

if size(x,1) == 1
   x = x(:);
   y = y(:);
   upper = upper(:);
   lower = lower(:);
end

N = size(x,2);
xy_line = cell(N,1);
xy_ci = cell(N,1);

clr_tr = 'rgba(255,255,255,0)';

for i = 1 : N

%     clr = round(colours(i,:) * 255);
%     clr_line = sprintf('rgb(%d,%d,%d)',clr(1),clr(2),clr(3));
%     clr_fill = sprintf('rgba(%d,%d,%d,0.5)',clr(1),clr(2),clr(3));
%     clr_tr = 'rgba(255,255,255,0)';

    clr_line = colours{i};
%     clr_fill = strrep(clr_line, '#', '#64');
    if clr_line(1) == '#'
        clr_fill = [clr_line dec2hex(alpha)];
    else
        clr_fill = strrep(strrep(clr_line,'rgb','rgba'),')',',0)');
    end

    xi = x(:,i);
    yi = y(:,i);
    upperi = upper(:,i);
    loweri = lower(:,i);
    
    xy_ci_i = struct(...
            'x', [xi;flip(xi)], ...
            'y', [upperi;flip(loweri)], ...
            'name', name{i}, ...
            'visible', 1, ...
            'fill', 'toself', ...
            'mode', 'lines', ...
            'fillcolor', clr_fill, ...
            'line', struct('color', clr_tr), ...
            'showlegend', false ...
            );

    xy_line_i = struct(...
            'x', xi, ...
            'y', yi, ...
            'name', name{i}, ...
            'visible', 1, ...
            'mode', 'lines+markers', ...
            'marker', struct(...
                            'symbol', 'circle', ...
                            'size', 13, ...
                            'color', clr_line ...
                            ), ...
            'line', struct('color', clr_line, ...
                           'width', 5), ...
            'showlegend', true ...
            );
     
     xy_line(i) = {xy_line_i};
     xy_ci(i) = {xy_ci_i};
        
end
    
end

