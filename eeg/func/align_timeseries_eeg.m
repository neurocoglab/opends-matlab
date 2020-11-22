function [ idx_map, d ] = align_timeseries_eeg( t1, t2 )
% Finds the closest corresponding time points in t2 for each time point in
% t1. 
%
% Arguments:
%  t1:  N1x1 vector
%  t2:  N2x1 vector
% 
% Returns:
%  idx_map:     Mx2 matrix with t1 indices and mapped t2 indices. Where
%                multiple points in t1 map to the same point in t2, only
%                the closest point is retained.
%  d:           The distances between the mapped points


    idx2 = nan(size(t1));
    dd = nan(size(t1));

    for i = 1 : length(t1)
        [dd(i),idx2(i)] = min(abs(t2-t1(i)));
    end

    ii = unique(idx2);
    idx_map = nan(length(ii),2);
    d = nan(size(ii));
    
    for i = 1 : length(ii)
        idxi = find(idx2==ii(i));
        [d(i),idxm] = min(dd(idxi));
        idx_map(i,1) = idxi(idxm);
        idx_map(i,2) = ii(i);
    end
    
end

