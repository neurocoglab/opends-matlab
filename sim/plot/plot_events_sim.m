function [ h ] = plot_events_sim( params, data, events, out2file )

    if nargin < 4
        out2file = false;
    end

    time_min = data.eye.t / 60000;
            
    gap_clr = [0.95 0.9 0.9];
    round_clr = [0.8 0.8 0.9];
    round_line_color = [1 0 0];
    overtake_color = [0 .5 0];
    noisy_colour = [0.9 0.9 0.9];
    left_change_color = [0.78 0 0.78];
    right_change_color = [0.82 0 0.3];
    pass_clr = [0.9 0.9 0.95];
    baseline_clr = [0.9 0.95 0.9];
    saccade_line_color = [0.3 0.3 0.3];
    saccade_rate_clr = [1 0 1];
    
    if out2file
        h = figure('visible','off');
    else
        h = figure;
    end
    set(h, 'Color', 'w');
    
    % Plot rounds background
%     if any(strcmp(events,'Rounds'))
%         
%         x1=0;
%         paint = false;
%         for j = 1 : size(results.sim2track.cycle_times,1)+1
%             
%             if j > size(results.sim2track.cycle_times,1)
%                 x2 = max(results.t);
%             else
%                 x2 = results.sim2track.cycle_times(j);
%             end
%             x2 = x2 / 60000;
% 
%             if paint
%                 rectangle('Position',[x1 -500 x2-x1 1000], ...
%                          'EdgeColor','w', ...
%                          'FaceColor', round_clr); 
%             end
%             
%             x1 = x2;
%             paint = ~paint;
%         end
%     end

    % Plot baseline periods
    if any(strcmp(events,'Baseline'))
        
       intervals = data.sim.sim2track.baseline;
       
       for j = 1 : size(intervals,1)
          x1 = intervals(j,1) / 60000;
          x2 = intervals(j,2) / 60000;
          
          if x2 > x1
          rectangle('Position',[x1 -500 x2-x1 1000], ...
                             'EdgeColor','w', ...
                             'FaceColor', baseline_clr); 
          end
       end
         
    end
    
    
    % Plot lane changes
    if any(strcmp(events,'LaneChanges'))
        
        max_interval = 30000; % Half a minute
        
        left = [data.sim.sim2track.left_change_times, ...
                true(length(data.sim.sim2track.left_change_times),1)];
        right = [data.sim.sim2track.right_change_times, ...
                 false(length(data.sim.sim2track.right_change_times),1)];
             
        left_right = [left;right];
        [~,idx] = sort(left_right(:,1));
        left_right = left_right(idx,:);
        
        this_left = -1;
        this_right = -1;
        for j = 1 : length(left_right)
           if left_right(j,2)
              % Is change to left 
              this_left = left_right(j,1);
           else
              % Is change to right
              this_right = left_right(j,1);
              if this_left > 0 && this_right-this_left < max_interval
                  % Valid passing segment, plot
                  x1 = this_left / 60000;
                  x2 = this_right / 60000;
                  if x2 > x1
                  rectangle('Position',[x1 -500 x2-x1 1000], ...
                             'EdgeColor','w', ...
                             'FaceColor', pass_clr); 
                  end
              end
              % Reset
              this_left = -1;
              this_right = -1;
           end
        end
        
%         for j = 1 : length(results.sim2track.left_change_times)
%             x1 = results.sim2track.left_change_times(j);
%             x1 = x1 / 60000;
%             line([x1 x1], [-500 500], ...
%                  'Color',left_change_color);
%             
%         end
%         
%         for j = 1 : length(results.sim2track.right_change_times)
%             x1 = results.sim2track.right_change_times(j);
%             x1 = x1 / 60000;
%             line([x1 x1], [-500 500], ...
%                  'Color',right_change_color);
%             
%         end
        
        
    end
    
    % Plot gaps
    for j = 1 : size(data.eye.tgap,1)
       x1 = time_min(data.eye.tgap(j,1));
       x2 = x1 + data.eye.tgap(j,2) / 60000;
       if x2 > x1
       rectangle('Position',[x1 -500 x2-x1 1000], ...
                 'EdgeColor','w', ...
                 'FaceColor', gap_clr); 
       end
    end


    hold on;
    
    % Plot diameter
    pdiam = data.eye.diam;
    if isfield(data.eye, 'blinks')
       pdiam = data.eye.blinks.diam; 
    end
    plot(time_min, pdiam, 'Color', noisy_colour);
    xsm = smooth(pdiam,250,'sgolay');
    plot(time_min, xsm, 'b');
    
    % Really smooth
    xxsm = smooth(pdiam,20000,'moving');
    plot(time_min, xxsm, 'Color', [.3 .3 .3], 'LineWidth', 1.5);
    
    minx = min(xsm); maxx = max(xsm);
    minx = minx-minx*0.1; maxx = maxx+maxx*0.1;
    
    ymax = max(xsm);
    
     % Plot saccade rate
    if any(strcmp(events,'SaccadeRate'))
        sr = zscore(data.eye.saccades.saccade_rate) + mean(pdiam);
        plot(data.eye.t / 60000, ...
             sr, ...
             'Color', saccade_rate_clr, ...
             'LineWidth', 1.5);
    end
    
    % Plot round lines
    if any(strcmp(events,'Rounds'))
        
        % First cycle
        hh = text(0.05, double(maxx-0.01*maxx),'1.1');
            set(hh,'Color', round_line_color);
            set(hh,'FontWeight', 'bold');
        
        % Cycles
        for j = 1 : length(data.sim.sim2track.cycle_times)
            x1 = data.sim.sim2track.cycle_times(j);
            x1 = x1 / 60000;
           
            line([x1 x1], [-500 500], ...
                 'Color', round_line_color, ...
                 'LineWidth', 1.5);
           
            hh = text((x1+0.05), double(maxx-0.01*maxx),sprintf('%d.1',j+1));
            set(hh,'Color', round_line_color);
            set(hh,'FontWeight', 'bold');
             
        end
        
        % Repeats
        cycle = 1; repeat = 2;
        for j = 1 : length(data.sim.sim2track.repeat_times)
           x1 = data.sim.sim2track.repeat_times(j);
           
           while cycle <= length(data.sim.sim2track.cycle_times) && ...
                 data.sim.sim2track.cycle_times(cycle) < x1
               cycle = cycle + 1;
               repeat = 2;
           end
           
           x1 = x1 / 60000;
           line([x1 x1], [-500 500], ...
                 'Color',round_line_color, ...
                 'LineStyle', '--', ...
                 'LineWidth', 1.5);
           
           hh = text((x1+0.05), double(maxx-0.01*maxx),sprintf('%d.%d',cycle,repeat));
           set(hh,'Color', round_line_color);
           set(hh,'FontWeight', 'bold');
           
           repeat = repeat + 1;
        end
        
        
    end
    
    % Plot overtake events
    if any(strcmp(events,'Overtakes'))
        
        for j = 1 : length(data.sim.sim2track.overtake_times)
            x1 = data.sim.sim2track.overtake_times(j);
            x1 = x1 / 60000;
           line([x1 x1], [-500 500], ...
                 'Color',overtake_color);
            
        end
        
    end
    
    
    % Plot saccades
    if any(strcmp(events,'Saccades'))
        for i = 1 : size(data.eye.saccades.saccades,1)
           i_m = data.eye.saccades.saccades(i,2);
           x1 = time_min(i_m);
           line([x1 x1], [-500 500], ...
                 'Color', saccade_line_color, ...
                 'LineWidth', 0.5);
            
        end
    end

%     grid on;
    
    hh = title('Pupil diameter with simulation events');
    set (hh, 'FontSize', 16);
    
    hh = xlabel('Time (min)');
    set (hh, 'FontSize', 14);
    hh = ylabel('Pupil diameter (mm)');
    set (hh, 'FontSize', 14);
    
    ylim([minx,maxx]);
    h.Position=([400 400 1500 500]);
%     resize_window(h, [1500 500], [400 400]);
    
    if out2file
        outdir = sprintf('%s/%s/figures', params.io.output_dir, data.subject);
        if ~exist(outdir,'dir')
           mkdir(outdir); 
        end
        saveas(h, sprintf('%s/events_sim.fig', outdir));
        xlim([0 3]);
        saveas(h, sprintf('%s/events_sim.png', outdir));
        close(h);
    end

end

