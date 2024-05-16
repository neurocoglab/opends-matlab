function [ cycle ] = get_cycle( sim2track, t )
%GET_CYCLE Gets the cycle for time t

cycle = nan;

idx = find(t >= [0;sim2track.cycle_times], 1, 'last');
if ~isempty(idx)
   cycle = idx;
end

end

