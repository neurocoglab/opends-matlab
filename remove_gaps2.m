function [ data ] = remove_gaps2( data, params )
%
% Remove gaps of continuous segments below a threshold
% as well as combining gaps separated by no more than a
% certain interval.
%

N = length(data.t);
isgap = zeros(N,1);

x = data.(params.criterion);

tstep = 1/data.Fs*1000;

i = 1;
while i < N
   
   if x(i) < params.gapthres
      
      i_start = i;
      i_test = i + params.gapmin-1;
      if i_test > N, i_test = N; end
      
      if sum(x(i:i_test)<params.gapthres) > (i_test - i) * 0.9 % at least 90% below threshold
         
         i = i_test + 1;
         buffer = 0;
         while (i < N && x(i) < params.gapthres) || buffer < 3
             i = i + 1;
             if (i < N && x(i) < params.gapthres)
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

dt = diff(data.t);
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
data.t = data.t(~isgap);
data.pos_left_x = data.pos_left_x(~isgap);
data.pos_left_y = data.pos_left_y(~isgap);
data.diam_left = data.diam_left(~isgap);
  

% dt = diff(data.t);
% idx = dt > tstep + 1;
% dt = dt(idx);
% tgap = find(idx);
% 
% isgap = zeros(length(x),1);
% for i = 1 : length(tgap)
%     idx = tgap(i) + 1;
%     while idx < N && x(idx) < params.gapthres
%         isgap(idx) = 1;
%         idx = idx + 1;
%     end
%     idx = tgap(i);
%     while idx > 0 && x(idx) < params.gapthres
%         isgap(idx) = 1;
%         idx = idx - 1;
%     end
% end
% 
% for i = 1 : length(tgap)
%     idx = tgap(i) + 1;
%     while idx < N && x(idx) > params.gapmax
%         isgap(idx) = 1;
%         idx = idx + 1;
%     end
%     idx = tgap(i);
%     while idx > 0 && x(idx) > params.gapmax
%         isgap(idx) = 1;
%         idx = idx - 1;
%     end
% end
% 
% data.t = data.t(~isgap);
% data.pos_left_x = data.pos_left_x(~isgap);
% data.pos_left_y = data.pos_left_y(~isgap);
% data.diam_left = data.diam_left(~isgap);

dt = diff(data.t);
idx = dt > tstep * params.gapmin;
dt = dt(idx);
data.tgap = [find(idx), dt];

end

