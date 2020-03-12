function [ times ] = interpolate_log_times_sim( M, sim_times, extrapolate )
%%%%%%%%%%%%%
% interpolate_log_times Converts simulation times to eye tracker times by 
% interpolating based on the matrix M, a row-sorted matrix mapping these times 
% based on simulation events. All simulation times to be mapped must lie 
% within the time points specified by M (i.e., the first and last rows).
%
% Note: This is necessary because timers for both systems are not
% synchronous (discrepancy of ~150 ms over 40 mins)
%

if nargin < 4
   extrapolate = false; 
end

% idx_ttime = find(strcmp(hdr, 'TrackerTime'));
% idx_stime = find(strcmp(hdr, 'SimTime'));

ttimes = M.TrackerTime;
stimes = M.SimTime;

% ttimes = double(cell2mat(M(:,idx_ttime)));
% stimes = cell2mat(M(:,idx_stime));

outcount = 0;
times = nan(length(sim_times),1);

for i = 1 : length(sim_times)
    ti = sim_times(i);
    idx = find(stimes > ti, 1);
    
    if (isempty(idx) || idx == 1) && ~extrapolate
       outcount = outcount + 1;
    else
       if isempty(idx), idx = length(stimes); end
       xs = stimes(idx-1);
       xe = stimes(idx);
       di = double(ti - xs) / double(xe - xs);
       
       ys = ttimes(idx-1);
       ye = ttimes(idx);
       
       times(i) = round(ys + (ye - ys) * di);
        
    end

end

if outcount > 5
   warning(['%d time points could not be extrapolated.' ...
            ' These entries will have nan values.'], outcount);
end

end

