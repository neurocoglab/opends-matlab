function [  ] = resize_window( h, size, origin )
%RESIZE_WINDOW Summary of this function goes here
%   Detailed explanation goes here

p = get(h, 'Position');
p(3:4)=size;
if exist('origin','var'), p(1:2) = origin; end;
set(h, 'Position', p);

end

