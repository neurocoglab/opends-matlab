function [ data ] = process_fixation_events( data, params )

%PROCESS_FIXATION_EVENTS Detects gaze intersections with AOIs
%
%   Areas of interest (AOIs) are defined by params.sim.events.fixations.aoi_file
%    This is a CSV file which specifies the name, round, distances, and rectangle
%    defining the AOI. Multiple distances/AOIs can be defined for the same
%    name. In this case, AOIs will be linearly interpreted based upon the
%    current position if the driver vehicle relative to the two distances
%    between which it is currently situated.
%
%   This function defines fixation window as the period of time during which the
%    gaze position is inside the AOI rectangle, subject to the parameters
%    outlined below.
%
%   Constraint/parameters are defined by params.eye.events.fixation:
%     1. min_duration: The minimal duration in seconds for which to accept
%                      a fixation
%     2. max_gap: The maximal proportion of the fixation window which is a
%                 gap (i.e., no signal detected)
%     3. smooth_span: The window (in samples) defining a moving average 
%                      smoothing of the gaze position. If zero, no smoothing 
%                      will be applied
%     4. max_errant: The maximal duration, in seconds, for which the gaze 
%                    position is errabnt (momentarily outside the AOI).
%     
%   Output is a table in data.events.fixations, which includes
%     1. Time of onset/offset
%     2. Track distance of onset/offset
%     3. Name of AOI
%     4. Round/repeat
%     

% 1. Load AOIs
aois_table = readtable(params.sim.events.fixations.aoi_file);

% 2. Smooth gaze position if necessary
eye_x = data.eye.pos_x;
eye_y = data.eye.pos_y;
if params.eye.events.fixations.smooth_factor > 0
   eye_x = smooth(eye_x, params.eye.events.fixations.smooth_span, 'moving');
   eye_y = smooth(eye_y, params.eye.events.fixations.smooth_span, 'moving'); 
end

% 3. TODO: get lane information; AOIs will be lane-dependent

idx_in = [];
gap_duration = 0;

% Create empty results table and append to it
varnames = {'Name', 'Round', 'Repeat', 'TimeStart', 'Duration', 'DistanceStart', 'Distance'};
results_table = cell2table(cell(0,7), 'VariableNames', varnames);

% 3. For each round
for r = 1 : N_rounds
%   Get AOIs for this round
    aois = aois_table(aois_table.round == r,:);
    for rp = 1 : N_repeat
%   Iterate through repeats/distances
        last_inside = 0;
        fixation_start = 0;
        duration_in = 0;
        duration_out = 0;
        for s = 1 : N_samples
            d = data.eye.distance(s);
            t = data.eye.t(s);
            if isempty(idx_in)
%               If gaze_position is in AOI
                for a = 1 : height(aois)
                    if d >= aois.dist_start(a) && d <= aois.dist_end(a)
                       idx_in = a;
                       window_start_t = t;
                       window_start_d = d;
                       window_name = aois.name{a};
                       aoi_start = [aois.rect_x(a) aois.rect_y(a) ...
                                    aois.rect_width(a) aois.rect_height(a)];
                       % Is there a next window?
                       if idx_in < height(aois) && strcmp(aois.name(a+1), window_name)
                           aoi_stop = [aois.rect_x(a+1) aois.rect_y(a+1) ...
                                        aois.rect_width(a+1) aois.rect_height(a+1)];
                       else
                           aoi_stop = [];
                       end
                       
                       continue;
                    end
                end
            end
            
            if ~isempty(idx_in)
                row = aois_table(idx_in,:);
                % If current_distance is in AOI time window:
                if ~(d >= row.dist_start && d <= row.dist_end)
                    % We are outside the time window
                    if isempty(aoi_stop)
                    % If next AOI is different, terminate fixation: 
                                        % if constraints are met, 
                                        % log fixation window
                        if fixation_start > 0
                            t_start =  data.eye.t(fixation_start);
                            d_start = data.eye.distance(fixation_start);
                            duration = t - t_start;
                            width = d - d_start;
                            if width > params.eye.events.fixations.min_duration && ...
                               gap_width / width < params.eye.events.fixations.max_gap

                                row_out = cell2table({name, r, rp, t_start, duration ...
                                                      d_start, width}, ...
                                                      'VariableNames', varnames);
                                results_table = [results_table;row_out];
                                
                            end
                        end
                        
                        last_inside = 0;
                        fixation_start = 0;
                        duration_in = 0;
                        duration_out = 0;
                        idx_in = [];
                        continue;
                    else
                        % Otherwise, we just need to change the AOI rectangle
                        aoi_start = aoi_stop;
                        % Is there a next window?
                        if idx_in < height(aois) && strcmp(aois.name(idx_in+1), window_name)
                           aoi_stop = [aois.rect_x(idx_in+1) aois.rect_y(idx_in+1) ...
                                        aois.rect_width(idx_in+1) aois.rect_height(idx_in+1)];
                        else
                           aoi_stop = [];
                        end
                        idx_in = idx_in+1;
                        row = aois_table(idx_in,:);
                    end
                end
                
                % Interpolate current AOI
                if ~isempty(aoi_stop)
                    % Need to interpolate the rectangle
                    dt = (d - window_start_d) / (row.dist_end / row.dist_start);
                    aoi_current = aoi_start + dt * (aoi_stop - aoi_start); 
                else
                    aoi_current = aoi_start;
                end
                    
                % Is the gaze inside the AOI?            
                if is_inside( aoi_current, gaze_x(s), gaze_y(s) )
                    % Is this the start of a fixation?
                    if fixation_start == 0
                        fixation_start = s;
                    end
                    
                    duration_in = duration_in + 1/data.eye.Fs;
                    last_inside = s;
                    duration_out = 0;
                elseif duration_in > 0
                    % if gaze_position is outside AOI, test whether to
                    % terminate
                    duration_out = duration_out + 1/data.eye.Fs;
                    if duration_out > params.eye.events.fixations.max_errant
                        if last_inside > 0
                            % We have a potentially valid fixation event, so add to log
                            t_start =  data.eye.t(fixation_start);
                            d_start = data.eye.distance(fixation_start);
                            duration = t - t_start;
                            width = d - d_start;
                            if width > params.eye.events.fixations.min_duration && ...
                               gap_width / width < params.eye.events.fixations.max_gap

                                row_out = cell2table({name, r, rp, t_start, duration ...
                                                      d_start, width}, ...
                                                      'VariableNames', varnames);
                                results_table = [results_table;row_out];

                            end
                        end
                        
                        last_inside = 0;
                        fixation_start = 0;
                        duration_in = 0;
                        duration_out = 0;
                        idx_in = [];
                        
                    end
                end

            end
        end
    end
end


end

