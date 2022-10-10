function [ clr ] = lighten_colour( clr, factor )
%LIGTHEN_COLOUR Lightens clr by interpolated towards white, by
%               the specified factor (0 to 1)

w = [1 1 1];
d = w - clr;
d = d * factor;
clr = clr + d;

end

