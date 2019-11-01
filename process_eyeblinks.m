%%
% Finds and removes eye blinks in eye tracking data
%
% Author: Andrew Reid
% Copyright (c) Andrew Reid, Univeristy of Nottingham, 2018
%
% Inputs:
%   data:     A struct containing the pupil data, with fields:
%                {variable-name}: A variable containing important eye data (e.g., 
%                                 gaze position, pupil diameter
%                t: Time (ms) relative to start
%                Fs: Sample frequency (Hz)
%                tgap: Represents gaps in the data; tgap(1) is the index of the gap;
%                      tgap(2) is the duration of the gap in ms. Note, function will
%                      process segments between gaps separately.
%
%   params:   A struct specifying parameters of the operation:
%                criterion: Name of the 'data' field to be used as the criterion for 
%                           eye-blink detection (typically this is pupil diameter)
%                smooth_width: Determines how much to smooth the data. See the 'smooth'
%                              function for details.
%                thres: Threshold for identifying a blink event; thres(1) identifies the
%                       negative-going onset, and thres(2) identifies the positive-going
%                       reopening.
%                interval: Defines the time interval over which to interpolate across the
%                          eyeblink events. interval(1) is time before onset and 
%                          interval(2) is time after reopening.
%                window: Time window for computing the blink rate as a moving average.
%                        Specified in ms.
%                maxblink: Maximal duration of a blink event; blinks longer than this will 
%                          be ignored.
% 
% Outputs:
%
%    results: A struct containing the results of the blink removal:
%                 blinks: Logical array where 1 indicates the midpoint of a blink event
%                 blink_ints: Actual blink intervals, where blink_ints(1) is the index of
%                             the blink, and blink_ints(2) is its duration in time points
%                 blink_rate: The blink rate computed as a moving average with time window
%                             specified by params.window
%                 d_x: The smoothed derivative of the criterion variable (i.e., velocity)
%                 intervals: The intervals on which blinks were detected
%


function [ results ] = process_eyeblinks (data, params )
       
x = data.(params.criterion);  % Usually, pupil diameter
N = length(x);

% Trim data
delta_t = 1 / data.Fs * 1000; % data.t(2)-data.t(1);

results.t = data.t; 
results.t_start = data.t_start;
results.tgap = data.tgap;

results.blinks = zeros(N,1);
results.blink_rate = zeros(N,1);
results.intervals = cell(0);
results.pct_fixed = cell(0);
results.pos_left_x = data.pos_left_x; % nan(N,1);
results.pos_left_y = data.pos_left_y; %nan(N,1);
results.diam_left = data.diam_left; %nan(N,1);
results.d_x = zeros(N,1);
results.blink_ints = cell(0);

% Process each segment seperately
idx = 1;
for j = 1 : size(data.tgap,1) + 1
   if j > size(data.tgap,1)
      if idx >= N, break; end
      idx2 = N;
   else
      idx2 = data.tgap(j,1);
   end

   if idx ~= idx2
       xx = x(idx:idx2);

       xx = smooth(xx,'moving',params.smooth_width);
       data_seg = [{data.pos_left_x(idx:idx2)}, ...
                   {data.pos_left_y(idx:idx2)}, ...
                   {data.diam_left(idx:idx2)}];

       [results.blinks(idx:idx2), results.blink_rate(idx:idx2), ...
        blink_ints, intervals_j, pct_fixed_j, d_x, ...
        data_seg_out] = ...
           process_segment(xx, data_seg, [idx idx2]);

       results.intervals(j) = {intervals_j};
       results.pct_fixed(j) = {pct_fixed_j};
       results.pos_left_x(idx:idx2) = data_seg_out{1};
       results.pos_left_y(idx:idx2) = data_seg_out{2};
       results.diam_left(idx:idx2) = data_seg_out{3};
       results.d_x(idx:idx2) = d_x;
       if ~isempty(blink_ints)
           results.blink_ints(end+1) = {blink_ints};
       end
   end
   
   idx = idx2 + 1; % data.eye.tgap(i,2);
end

% Interpolate over gaps
buffer = round(params.gapbuffer / data.Fs * 1000 / 2);

for j = 1 : size(data.tgap,1)
    try
        idx = data.tgap(j,1);
        idx0 = max(1,idx-buffer);
        idx1 = min(N,idx+buffer);
        results.diam_left(idx0:idx1) = ...
            linterp([idx0 idx1],[results.diam_left(idx0) results.diam_left(idx1)],idx0:idx1);
    catch err
        fprintf('Error: %s', err)
    end
end

% patch up nans; TODO: why are these nans?
idx = find(isnan(results.pos_left_x));
results.pos_left_x(idx) = results.pos_left_x(idx-1);
results.pos_left_y(idx) = results.pos_left_y(idx-1);
results.diam_left(idx) = results.diam_left(idx-1);

%TODO: Combine blink ints into one list

    function [ blinks, blink_rate, blink_ints, intervals, pct_fixed, d_x, data_seg_out ] = ...
                                                process_segment( seg, data_seg, seg_idx )
        d_x = [0;diff(seg) / delta_t];
        %results.diff = d_x;
        data_seg_out = data_seg;

        window = ceil(params.window / delta_t);
        if mod(window,2) == 0, window = window+1; end

        Ns = length(seg);
        intervals = seg < params.absthres;
        
        int_left = ceil(params.interval(1) / delta_t);
        int_right = ceil(params.interval(2) / delta_t);

        % Identify blink intervals
        % Overlapping intervals will be merged
        blinks = zeros(Ns,1);
        blink_ints = [0 0];
        i = 1;
        jj = 1;
        while i <= Ns

%             try
            if d_x(i) < params.thres(1)
                i_b0 = i;
                i_s = i - int_left; if i_s < 1, i_s = 1; end
                i = i + 1;

                ismax = false;
                while i <= Ns && d_x(i) < params.thres(2) && ~ismax
                    i = i + 1;
                    ismax = (i - i_b0) > params.maxblink;
                end

				if ~ismax
					i = i + int_right;
					if i > Ns, i = Ns; end
					i_b1 = i;
					i_e = i;  %+ int_right; if i_e > Ns, i_e = Ns; end 

					intervals(i_s:i_e) = true;

					blinks(round((i_b1+i_b0)/2))=1;
					blink_ints(jj,:) = [i_s+seg_idx(1) i_e-i_s-1];
					jj = jj + 1;
				end
            end
%             catch err
%                a = 0; 
%             end
            i = i + 1;
        end
        
        % Get blinks as intervals
        if jj < 2
            blink_ints = [];
        end

        % Compute weighted moving average blink rate
        blink_rate = zeros(Ns,1);
        w = normpdf(-4:8/(window-1):4,0,1);

        for i = 1 : Ns

            i0=1; i1=window;
            i_s = i - (window-1)/2;
            if i_s < 1, i0 = 2-i_s; i_s = 1; end
            i_e = i + (window-1)/2;
            if i_e > Ns, i1 = window-(i_e-Ns); i_e = Ns; end

            wi = w(i0:i1)';
            blink_rate(i) = sum(blinks(i_s:i_e).*wi) / sum(wi);

        end

        pct_fixed = sum(intervals)/Ns;

        % For all time series, interpolate across intervals
        i = 1;
        while i <= Ns
           if intervals(i)
               i_s = i;
               i = i + 1;
               while i <= Ns && intervals(i) 
                   i = i + 1;
               end
               if i == Ns
                   i_e = i;
               else
                   i_e = i - 1;
               end

               % Interpolate
               for s = 1 : 3
                   ds = data_seg_out{s};
                   yy = linterp([i_s i_e],[ds(i_s) ds(i_e)], i_s:i_e);
                   
                   ds(i_s : i_e) = yy;
%                    if ~isempty(params.smooth)
%                       ds = smooth(ds,params.smooth,params.smooth_width); 
%                    end
                   data_seg_out(s) = {ds};
               end

           end
           
           i = i + 1;
        end
        
       % Smooth all segments
       for itr = 1 : 1
           for s = 1 : 3
               ds = data_seg_out{s};
               if ~isempty(params.smooth)
                  ds = smooth(ds,'moving',params.smooth_width); 
               end
               data_seg_out(s) = {ds};
           end
       end
        
    end

end

