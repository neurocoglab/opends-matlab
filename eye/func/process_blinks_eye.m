%%
% Finds and removes eye blinks in eye tracking data.eye
%
% Author: Andrew Reid
% Copyright (c) Andrew Reid, Tilburg University, 2024
%
% Inputs:
%   data.eye:     A struct containing the pupil data.eye, with fields:
%                {variable-name}: A variable containing important eye data.eye (e.g., 
%                                 gaze position, pupil diameter
%                t: Time (ms) relative to start
%                Fs: Sample frequency (Hz)
%                tgap: Represents gaps in the data.eye; tgap(1) is the index of the gap;
%                      tgap(2) is the duration of the gap in ms. Note, function will
%                      process segments between gaps separately.
%
%   params.eye.blinks:   A struct specifying parameters of the operation:
%                criterion: Name of the 'data.eye' field to be used as the criterion for 
%                           eye-blink detection (typically this is pupil diameter)
%                smooth_width: Determines how much to smooth the data.eye.blinks. See the 'smooth'
%                              function for details.
%                thres: Threshold for identifying a blink event; thres(1) identifies the
%                       negative-going onset, and thres(2) identifies the positive-going
%                       reopening.
%                interval: Defines the time interval over which to interpolate across the
%                          eyeblink events. interval(1) is time before onset and 
%                          interval(2) is time after reopening. Specified in ms.
%                window: Time window for computing the blink rate as a moving average.
%                        Specified in ms.
%                maxblink: Maximal duration of a blink event; blinks longer than this will 
%                          be ignored. Specified in ms.
% 
% Outputs:
%
%    data: A updated data struct (new field "blinks") containing the results of the blink removal:
%                 blinks: Logical array where 1 indicates the midpoint of a blink event
%                 blink_ints: Actual blink intervals, where blink_ints(1) is the index of
%                             the blink, and blink_ints(2) is its duration in time points
%                 blink_rate: The blink rate computed as a moving average with time window
%                             specified by params.eye.blinks.window
%                 d_x: The smoothed derivative of the criterion variable (i.e., velocity)
%                 intervals: The intervals on which blinks were detected
%                 diam, 
%


function [ data ] = process_blinks_eye ( data, params )
       
x = data.eye.(params.eye.blinks.criterion);  % Usually, pupil diameter
x(isnan(x)) = 0;
N = length(x);

delta_t = 1 / data.eye.Fs * 1000;
max_blink = params.eye.blinks.maxblink * delta_t;

data.eye.blinks = [];
data.eye.blinks.blinks = zeros(N,1);
data.eye.blinks.blink_rate = zeros(N,1);
% data.eye.blinks.intervals = cell(0);
data.eye.blinks.intervals = zeros(N,1);
data.eye.blinks.pct_fixed = cell(0);
data.eye.blinks.d_x = zeros(N,1);
data.eye.blinks.blink_ints = cell(0);
data.eye.blinks.params = params.eye.blinks;

% Process each segment seperately
idx = 1;
for j = 1 : size(data.eye.tgap,1) + 1
   if j > size(data.eye.tgap,1)
      if idx > N, break; end
      idx2 = N;
   else
      idx2 = data.eye.tgap(j,1);
   end

   if idx ~= idx2
       xx = x(idx:idx2);

       if params.eye.blinks.smooth
            xx = smooth(xx,'moving',params.eye.blinks.smooth_width);
       end
       data_seg = [{data.eye.pos_x(idx:idx2)}, ...
                   {data.eye.pos_y(idx:idx2)}, ...
                   {data.eye.diam(idx:idx2)}];

       [data.eye.blinks.blinks(idx:idx2), ...
        blink_ints, intervals_j, pct_fixed_j, d_x, ...
        data_seg_out] = ...
           process_segment(xx, data_seg, [idx idx2]);

       data.eye.blinks.intervals(idx:idx2) = intervals_j;
       data.eye.blinks.pct_fixed(j) = {pct_fixed_j};
       data.eye.blinks.pos_x(idx:idx2) = data_seg_out{1};
       data.eye.blinks.pos_y(idx:idx2) = data_seg_out{2};
       data.eye.blinks.diam(idx:idx2) = data_seg_out{3};
       data.eye.blinks.d_x(idx:idx2) = d_x;
       if ~isempty(blink_ints)
           data.eye.blinks.blink_ints(end+1) = {blink_ints};
       end
   else
       % end case
       data.eye.blinks.pos_x(idx) = data.eye.pos_x(idx);
       data.eye.blinks.pos_y(idx) = data.eye.pos_y(idx);
       data.eye.blinks.diam(idx) = x(idx);
   end
   
   idx = idx2 + 1;
end

data = get_blink_intervals( data );

% Add from gaps?
% if ~isempty(params.eye.blinks.from_gaps.width_lims)
%     
%     W = data.eye.tgap(:,2);
%     idx = find(W > params.eye.blinks.from_gaps.width_lims(1) & ...
%                W < params.eye.blinks.from_gaps.width_lims(2));
%     if ~isempty(idx)
%         new_intervals = zeros(length(idx),2);
%         new_blinks = zeros(length(idx),1);
%         for j = 1 : length(idx)
%             new_intervals(j,:) = [data.eye.tgap(1) data.eye.tgap(1)+1];
%             new_blinks(j) = data.eye.tgap(1);
%         end
%     end
% end

% Interpolate over gaps
if params.eye.blinks.gapbuffer > 0
    
    buffer = round(params.eye.blinks.gapbuffer / data.eye.Fs * 1000 / 2);

    for j = 1 : size(data.eye.tgap,1)
        try
            idx = data.eye.tgap(j,1);
            idx0 = max(1,idx-buffer);
            idx1 = min(N,idx+1+buffer);
            data.eye.blinks.diam(idx0:idx1) = ...
                linterp([idx0 idx1],[data.eye.blinks.diam(idx0) data.eye.blinks.diam(idx1)],idx0:idx1);
        catch err
            fprintf('Error: %s', err.message)
        end
    end
    
end

% Get blink rate
data = get_blink_rate( data );

% patch up nans if necessary
% idx = find(isnan(data.eye.blinks.pos_x));
% if ~isempty(idx)
%     data.eye.blinks.pos_x(idx) = data.eye.pos_x(idx-1);
%     data.eye.blinks.pos_y(idx) = data.eye.pos_y(idx-1);
%     data.eye.blinks.diam(idx) = data.eye.diam(idx-1);
% end

    function [ blinks, blink_ints, intervals, pct_fixed, d_x, data_seg_out ] = ...
                                                process_segment( seg, data_seg, seg_idx )
        d_x = [0;diff(seg) / delta_t];
        data_seg_out = data_seg;

        Ns = length(seg);
        intervals = false(Ns,1); % seg < params.eye.blinks.absthres;
        
        int_left = ceil(params.eye.blinks.interval(1) / delta_t);
        int_right = ceil(params.eye.blinks.interval(2) / delta_t);

        % Identify blink intervals
        % Overlapping intervals will be merged
        blinks = zeros(Ns,1);
        blink_ints = [0 0];
        i = 1;
        jj = 1;
        while i <= Ns

            if d_x(i) < params.eye.blinks.thres(1)
                i_b0 = i;
                i_s = i - int_left; if i_s < 1, i_s = 1; end
                i = i + 1;

                ismax = false;
                while i <= Ns && d_x(i) < params.eye.blinks.thres(2) && ~ismax
                    i = i + 1;
                    ismax = (i - i_b0) > max_blink;
                end

				if ~ismax
					i2 = i + int_right;
					if i2 > Ns, i2 = Ns; end
					i_b1 = i2;
					i_e = i2;

					intervals(i_s:i_e) = true;

					blinks(round((i_b1+i_b0)/2))=1;
					blink_ints(jj,:) = [i_s+seg_idx(1) i_e-i_s-1];
					jj = jj + 1;
				end
            end
            
            i = i + 1;
        end
        
        % Get blinks as intervals
        if jj < 2
            blink_ints = [];
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
                   data_seg_out(s) = {ds};
               end

           end
           
           i = i + 1;
        end
        
       % Smooth all segments
       for itr = 1 : 1
           for s = 1 : 3
               ds = data_seg_out{s};
               if ~isempty(params.eye.blinks.smooth)
                  ds = smooth(ds,'moving',params.eye.blinks.smooth_width); 
               end
               data_seg_out(s) = {ds};
           end
       end
        
    end

    function data = get_blink_rate( data ) 
        
        % Compute weighted moving average blink rate
        window = ceil(params.eye.blinks.window / delta_t);
        if mod(window,2) == 0, window = window+1; end
        blinks = data.eye.blinks.blinks;
        
        Ns = length(blinks);
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

        data.eye.blinks.blink_rate = blink_rate;
        
    end

    function data = get_blink_intervals( data )

        intervals = data.eye.blinks.intervals;
        blinks = zeros(length(intervals),1);
        dints = diff(double(intervals));
        ints = zeros(sum(dints==1),2);
        ii = 1; jj = 1; is_int = false;

        while ii <= length(dints)
           if is_int
               if dints(ii) == -1
                   ints(jj,:) = [i_s ii];
                   blinks(i_s+round((ii-i_s)/2)) = 1;
                   jj = jj + 1;
                   is_int = false;
                   i_s = -1;
               end
           else
               if dints(ii) == 1
                   i_s = ii; 
                   is_int = true;
               else
                   i_s = -1; 
                   is_int = false;
               end
           end
           ii = ii + 1;
           if ii > length(dints) && is_int
               ints(jj,:) = [i_s, length(intervals)];
               blinks(i_s+round((length(intervals)-i_s)/2)) = 1;
           end
        end
        
        data.eye.blinks.blinks = blinks;
        data.eye.blinks.intervals = ints;

    end

end

