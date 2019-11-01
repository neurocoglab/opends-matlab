% Init

Fs = 500;           % in hz
ts = LDiaXpx_sm2;
t = t_ts;
thres_low = -10;
thres_high = 10;
buffer = 10;        % in ms
method= 'spline'; % 'linear'
max_width=500;

bb = buffer / (1000/Fs);

% Smooth data

% hann=hanning(10);
% ts = conv(ts,hann,'same');

delta_t = t(2)-t(1);
d_ts = diff(ts) / delta_t * 100;
d_ts = [0 d_ts'];
d_ts(1:5)=0;
d_ts(length(d_ts)-5:end)=0;
f_ts = ts;

% Find blink timepoints
i = 1;
while i <= length(d_ts)

% Scan across time series derivative
% If d_ts(i) < thres_low, t2 = i, start
if d_ts(i) < thres_low
    t2 = i - buffer; if t2<1, t2=1; end
    
    isoffset = false;
    found = false;
    start = i;
    while (i-start < max_width) && i <= length(d_ts) && ~found
        if isoffset
           if d_ts(i) < thres_high
               found = true;
               t3 = i + buffer; if t3 > length(d_ts), t3=length(d_ts); end
           end
        elseif d_ts(i) > thres_high
           isoffset = true; 
        end
        
        if ~found
            i = i + 1;
        end
    end
    
    if found
        % If d_ts(j) > thres_high, t3 = j
        
        % t1 = t2-t3+t2; t4=t3-t2+t3
        t1 = t2-t3+t2; if t1<1, t1=1; end
        t4 = t3-t2+t3; if t4>length(d_ts), t4=length(d_ts); end
        
        % fit cubic spline
        if strcmp(method,'spline')
            yy = spline([t(t1) t(t2) t(t3) t(t4)],[f_ts(t1) f_ts(t2) f_ts(t3) f_ts(t4)],t(t1:t4));
            f_ts(t1:t4) = yy;
        elseif strcmp(method,'poly')
            p = polyfit([t(t1) t(t2) t(t3) t(t4)],[f_ts(t1) f_ts(t2) f_ts(t3) f_ts(t4)],2);
            yy = polyval(p,t(t1:t4));
            f_ts(t1:t4) = yy;
        else
            yy = linterp([t(t2) t(t3)],[f_ts(t2) f_ts(t3)], t(t2:t3));
            f_ts(t2:t3) = yy;
        end
        
    end
end

i = i + 1;
end
