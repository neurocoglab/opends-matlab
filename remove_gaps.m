function [ data ] = remove_gaps( data, params )
%
% Remove gaps of continuous segments below a threshold
% as well as combining gaps separated by no more than a
% certain interval.
%

N = length(data.t);
gaps = [];
isgap = zeros(N,1);

delta_t = 1/data.Fs * 1000;

x = data.(params.criterion);

i = 1;
while i < N
   
   if x(i) < params.gapthres
      
      i_start = i;
      i_test = i + params.gapmin-1;
      if i > 57300
         a=0; 
      end
      if i_test > N, i_test = N; end
      
      if sum(x(i:i_test)<params.gapthres) == i_test - i + 1
         
         i = i_test + 1;
         while i < N && x(i) < params.gapthres
             i = i + 1;
         end
         if i > N+1, i = N+1; end
         
         isgap(i_start : i - 1) = 1;
         gaps(end+1,:) = [i_start (i - i_start + 1)];
         
      end
      
   end
       
  i = i + 1;
end

   % Adjust gaps for removed items
   offsets = zeros(size(gaps,1),2);
   offset = 0;
   for i = 1 : size(gaps,1)
       offsets(i,1) = gaps(i,1);
       offsets(i,2) = offset;
       offset = offset + gaps(i,2);
   end

   % Sort gaps
   gaps(:,2) = gaps(:,2) * delta_t;
   data.tgap = [data.tgap;gaps];
   g = data.tgap(:,1);
   [~,idx] = sort(g);
   data.tgap = data.tgap(idx,:);
  
   % Adjust times
   j = 1;
   for i = 1 : size(data.tgap,1)
       while data.tgap(i,1) > offsets(j,1) && j < size(offsets,1)
           j = j + 1; 
       end
       if j <= size(offsets,1)
           data.tgap(i,1) = data.tgap(i,1) - offsets(j,2);
       end
   end
   
   
   % Remove gaps from time series
   data.t = data.t(~isgap);
   data.pos_left_x = data.pos_left_x(~isgap);
   data.pos_left_y = data.pos_left_y(~isgap);
   data.diam_left = data.diam_left(~isgap);
   
    
   % Merge gaps
   newgaps = data.tgap(1,:);
   offsets = [];
   offset = 0;
   o=1;
   isgap=zeros(length(data.t),1);
   for i = 2 : size(data.tgap,1)
       last_idx = data.tgap(i-1,1)+1;
       this_idx = data.tgap(i,1);
       if this_idx - last_idx < params.gapmerge
           offset = offset + this_idx - last_idx;
           offsets(o,:) = [this_idx, offset];
           isgap(last_idx:this_idx)=1;
           dt = data.t(this_idx)-data.t(last_idx);
       end
   end
   
   data.tgap




end

