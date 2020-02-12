function [ results ] = remove_gaps( data, params )
%
% Remove gaps of continuous segments below a threshold
% as well as combining gaps separated by no more than a
% certain interval.
%
%
% Inputs:
%   data:     A struct containing the pupil data, with fields:
%                {criterion variable}: variable containing gap criterion
%                                      information (typically this is pupil diameter)
%                t: Time (ms) relative to start
%                Fs: Sample frequency (Hz)
%
%   params:   A struct specifying parameters of the operation:
%                criterion:     Variable name in data struct to act as the
%                               criterion for determining gaps (typically this
%                               the pupil diameter)
%                variables:     Cell array of variable names from which to
%                               remove gaps; gaps will also be removed from
%                               data.t
%                gapthres:      Maximal threshold (in units of the criterion
%                               variable) at which to accept a gap (at
%                               least 90% of the interval must be below
%                               this value)
%                gapmin:        Miminal width (in samples) of a gap
%                gapmerge:      Maximal number of suprathreshold samples
%                               between two gaps to determine whether these
%                               gaps should be merged
%                gapthresmax:   If this field exists and is not empty, a
%                               values above this value will also be
%                               treated as gaps
% 
% 
% 
% Outputs:
%
%    results: A struct that is a copy of "data" with gaps removed from each
%             variable in params.variables. A new member called "tgap" will also be
%             added. This is an Nx2 matrix, where N=number of removed gaps, the
%             first row being the sample index of a gap and the second being its
%             width (in samples).
%                 
%

gapthresmax = Inf;
if isfield(params, 'gapthresmax') && ~isempty(params.gapthresmax)
   gapthresmax = params.gapthresmax;
end

results = data;

N = length(results.t);
isgap = zeros(N,1);

x = results.(params.criterion);

tstep = 1/results.Fs*1000;

i = 1;
while i < N
   
   if x(i) < params.gapthres || x(i) > gapthresmax
      
      i_start = i;
      i_test = i + params.gapmin-1;
      if i_test > N, i_test = N; end
      
      if sum(x(i:i_test) < params.gapthres | x(i:i_test) > gapthresmax) > (i_test - i) * 0.9 % at least 90% below threshold
         
         i = i_test + 1;
         buffer = 0;
         while (i < N && (x(i) < params.gapthres || x(i) > gapthresmax)) || buffer < 3
             i = i + 1;
             if (i < N && (x(i) < params.gapthres || x(i) > gapthresmax))
                 buffer = 0;
             else
                 buffer = buffer + 1;
             end
         end
         if i > N+1, i = N+1; end
         
         i_start = max(1,i_start-buffer);
%          i = min(N,i+buffer);
         
         isgap(i_start : i) = 1;
         
      end
      
   end
       
  i = i + 1;
end

% Merge gaps

dt = diff(results.t);
% if dt(1)>0, dt(1) = tstep * params.gapmin + 1; end
tgaps = dt > tstep * params.gapmin;
tgaps = [tgaps;0];

i = 1;
while i <= N
   
    i_restart = 0;
    if isgap(i) || tgaps(i)
       
        while i <= N && isgap(i), i = i + 1; end
        
        if i > N, break; end
        
        if tgaps(i), i = i + 1; end
  
        i_start = i;
        count=0;
        while i <= N && ~(isgap(i) || tgaps(i))
            count = count + 1;
            i = i + 1; 
        end
        
        if i <= N
            if tgaps(i)
                i_restart = i;
            end

            if count < params.gapmerge
               isgap(i_start:i) = true;
               i_restart = i;
            end
        end
        
    end
    
    if i_restart > 0
        i = i_restart;
    else
        i = i + 1;
    end
end

x = x(~isgap);

% Remove gaps from time series
results.t = results.t(~isgap);

for i = 1 : length(params.variables)
    variable = params.variables{i};
    V = results.(variable);
    results.(variable) = V(~isgap);
end

% results.pos_left_x = results.pos_left_x(~isgap);
% results.pos_left_y = results.pos_left_y(~isgap);
% V = results.(params.criterion);
% results.(params.criterion) = V(~isgap);

dt = diff(results.t);
idx = dt > tstep * params.gapmin;
dt = dt(idx);
results.tgap = [find(idx), dt];

end

