function [ Y ] = remove_eyeblinks (x, Fs, thres_min, thres_max, buffer, interval )
       
if nargin < 6
   interval = [0 0]; 
end
buffer = ceil(buffer * (Fs / 1000));
N = length(x);
Y = x;

% Stroll through time points
% If x(t) == 0, start blink

i = 1;
while i <= N

    if x(i) < thres_min
        i_s = i - interval(1); if i_s < 1, i_s = 1; end
        i = i + 1;
        
        b = buffer; b_start=-1;
        while i <= N && x(i)< thres_max && (x(i)<thres_min || b > 0)
            if x(i)>thres_min
                b=b-1;
                if b_start<0, b_start=i; end
            else
                b=buffer; b_start=-1;
            end
            i = i + 1;
        end
        
        if b_start > 0, i = b_start+1; end
        
        i = i - 1;
        i_e = i + interval(2); if i_e > N, i_e = N; end 
        yy = linterp([i_s i_e],[x(i_s) x(i_e)], i_s:i_e);
        Y(i_s:i_e) = yy;
        
    end
    i = i + 1;
end


end