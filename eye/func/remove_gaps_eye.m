function [ data ] = remove_gaps_eye( data, params )
%
% Remove gaps of continuous segments below a threshold
% as well as combining gaps separated by no more than a
% certain interval.
%

N = length(data.eye.t);
isgap = zeros(N,1);

x = data.eye.(params.eye.gaps.criterion);

tstep = 1/data.eye.Fs*1000;

i = 1;
while i < N
   
   if x(i) < params.eye.gaps.gapthres
      
      i_start = i;
      i_test = i + params.eye.gaps.gapmin-1;
      if i_test > N, i_test = N; end
      
      if sum(x(i:i_test)<params.eye.gaps.gapthres) > (i_test - i) * 0.9 % at least 90% below threshold
         
         i = i_test + 1;
         buffer = 0;
         while (i < N && x(i) < params.eye.gaps.gapthres) || buffer < 3
             i = i + 1;
             if (i < N && x(i) < params.eye.gaps.gapthres)
                 buffer = 0;
             else
                 buffer = buffer + 1;
             end
         end
         if i > N+1, i = N+1; end
         
         i_start = max(1,i_start-buffer);
         i = min(N,i+buffer);
         
         isgap(i_start : i) = 1;
         
      end
      
   end
       
  i = i + 1;
end

% Merge gaps

dt = diff(data.eye.t);
tgaps = dt > tstep * params.eye.gaps.gapmin;
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

            if count < params.eye.gaps.gapmerge
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

% Apply buffer in samples (specified in ms)
if params.eye.gaps.buffer > 0
    buffer = ceil(params.eye.gaps.buffer / (1000 / data.eye.Fs));
    dg = diff(isgap);
    inext = 0;
    for i = 1 : length(dg)
        if i > inext
            if dg(i) < 0
                idx_b = min(length(isgap),i+buffer);
                isgap(i:idx_b) = true;
                inext = idx_b;
            elseif dg(i) > 0
                idx_b = max(1,i-buffer);
                isgap(idx_b:i) = true;
            end
        end
    end
end

x = x(~isgap);

% Remove gaps from time series
data.eye.t = data.eye.t(~isgap);
data.eye.pos_x = data.eye.pos_x(~isgap);
data.eye.pos_y = data.eye.pos_y(~isgap);
data.eye.diam = data.eye.diam(~isgap);
  

% dt = diff(data.eye.t);
% idx = dt > tstep + 1;
% dt = dt(idx);
% tgap = find(idx);
% 
% isgap = zeros(length(x),1);
% for i = 1 : length(tgap)
%     idx = tgap(i) + 1;
%     while idx < N && x(idx) < params.eye.gaps.gapthres
%         isgap(idx) = 1;
%         idx = idx + 1;
%     end
%     idx = tgap(i);
%     while idx > 0 && x(idx) < params.eye.gaps.gapthres
%         isgap(idx) = 1;
%         idx = idx - 1;
%     end
% end
% 
% for i = 1 : length(tgap)
%     idx = tgap(i) + 1;
%     while idx < N && x(idx) > params.eye.gaps.gapmax
%         isgap(idx) = 1;
%         idx = idx + 1;
%     end
%     idx = tgap(i);
%     while idx > 0 && x(idx) > params.eye.gaps.gapmax
%         isgap(idx) = 1;
%         idx = idx - 1;
%     end
% end
% 
% data.eye.t = data.eye.t(~isgap);
% data.eye.pos_x = data.eye.pos_x(~isgap);
% data.eye.pos_y = data.eye.pos_y(~isgap);
% data.eye.diam = data.eye.diam(~isgap);

dt = diff(data.eye.t);
idx = dt > tstep * params.eye.gaps.gapmin;
dt = dt(idx);
data.eye.tgap = [find(idx), dt];

end

