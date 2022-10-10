function [ h ] = plot_ci_filled( x, y, upper, lower, colours, lighten, ...
                                 marker, marker_size, line_width )
%PLOT_CI_FILLED Plot a line with a filled confidence interval 
%  

if nargin < 9
   line_width = 2; 
end

if nargin < 8
   marker_size = 30; 
end

if nargin < 7
   marker = []; 
end

if nargin < 6
   lighten = 0.8; 
end

if nargin < 5
   colours = jet(size(x,1)); 
end

if sum(size(x)) == max(size(x))+1
    N_lines = 1;
    x = x(:);
    y = y(:);
    upper = upper(:);
    lower = lower(:);
    colours = colours(:)';
else
    N_lines = size(x,2);
end

% Plot bounds
for i = 1 : N_lines
    colour = colours(i,:);
    fclr = lighten_colour(colour,lighten);
    xx = x(:,i);
    uu = upper(:,i);
    ll = lower(:,i);
    hh = fill([xx;xx(end:-1:1)],[uu;ll(end:-1:1)],fclr);
    set(hh,'EdgeColor','none');
    alpha(hh,0.5);
    hold on;
end

% Plot lines
for i = 1 : N_lines
    colour = colours(i,:);
    xx = x(:,i);
    yy = y(:,i);
    if isempty(marker)
        h = plot(xx,yy);
    else
        h = plot(xx,yy, marker);
    end
    h.MarkerSize = marker_size;
    set(h,'Color',colour);
    set(h,'LineWidth',line_width);
    hold on;
end

 hold off;

end

