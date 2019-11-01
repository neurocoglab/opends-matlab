function [ h ] = plot_events2( results, params, events, eeg_data_in, channels, plot_title )

    time_sec = minutes(results.t / 60000);
            
    gap_clr = [0.95 0.9 0.9];
    round_clr = params.round_clr; 
    round_line_color = params.round_line_clr;
    overtake_color = params.overtake_clr;
    noisy_colour = [0.9 0.9 0.9];
    left_change_color = [0.78 0 0.78];
    right_change_color = [0.82 0 0.3];
    pass_clr = params.lane_changes_clr;
    baseline_clr = params.baseline_clr; % [0.9 0.95 0.9];
    saccade_line_color = params.saccade_line_clr;
    saccade_rate_clr = params.saccade_rate_clr;
    diam_clr = params.diameter_clr;
    
    t_end = minutes(results.sim2track.simended_time / 60000);
    idx_keep = time_sec <= t_end;
    time_sec = time_sec(idx_keep);
    
    if ~isempty(params.to_file)
        h = figure('visible', 'off');
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
    if params.plot_baseline && any(strcmp(events,'Baseline'))
        
       intervals = minutes(results.sim2track.baseline / 60000);
       
       for j = 1 : size(intervals,1)
          x1 = intervals(j,1); % / 60000;
          x2 = intervals(j,2); % / 60000;
          
          if x2 > x1
%           rectangle('Position',[x1 -500 x2-x1 1000], ...
%                              'EdgeColor','w', ...
%                              'FaceColor', baseline_clr); 
          fill([x1 x2 x2 x1], [-500 -500 1000 1000], baseline_clr, 'EdgeColor','w');
          hold on;               
          end
       end
         
    end
    
    
    % Plot lane changes
    if params.plot_lane_changes && any(strcmp(events,'LaneChanges'))
        
        max_interval = 30000; % Half a minute
        
        left = [results.sim2track.left_change_times, ...
                true(length(results.sim2track.left_change_times),1)];
        right = [results.sim2track.right_change_times, ...
                 false(length(results.sim2track.right_change_times),1)];
             
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
                  x1 = minutes(this_left / 60000);
                  x2 = minutes(this_right / 60000);
                  if x2 > x1
                      fill([x1 x2 x2 x1], [-500 -500 1000 1000], pass_clr, 'EdgeColor','w');
                      hold on;   
%                   rectangle('Position',[x1 -500 x2-x1 1000], ...
%                              'EdgeColor','w', ...
%                              'FaceColor', pass_clr); 
                  end
              end
              % Reset
              this_left = -1;
              this_right = -1;
           end
        end
        
        
        
    end
    
    % Plot gaps
    if params.plot_gaps
        for j = 1 : size(results.tgap,1)
           x1 = results.t(results.tgap(j,1));
           x2 = x1 + results.t(results.tgap(j,2));
           if x2 < results.sim2track.simended_time
               x1 = minutes(x1/60000);
               x2 = minutes(x2/60000);
               if x2 > x1
    %            rectangle('Position',[x1 -500 x2-x1 1000], ...
    %                      'EdgeColor','w', ...
    %                      'FaceColor', gap_clr); 
                    fill([x1 x2 x2 x1], [-500 -500 1000 1000], pass_clr, 'EdgeColor','w');
                    hold on;
               end
           end
        end
    end


    hold on;
    
    if params.plot_diameter
        
        diam = results.diam_left(idx_keep);
        if params.plot_diam_zscore
           diam = zscore(diam); 
        end
        
        % Plot diameter
        plot(time_sec, diam, 'Color', noisy_colour);
        xsm = smooth(diam,250,'sgolay');
        plot(time_sec, xsm, 'Color', diam_clr, 'LineWidth', params.diameter_linewidth);
        minx = min(xsm); maxx = max(xsm);
        minx = minx-abs(minx)*0.1; maxx = maxx+abs(maxx)*0.1;
    end
    
    % Really smooth
%     xxsm = smooth(results.diam_left,20000,'moving');
%     plot(time_sec, xxsm, 'Color', [.3 .3 .3], 'LineWidth', 1.5);
    
    
    
%     ymax = max(xsm);
    
     % Plot saccade rate
    if params.plot_saccade_rate
        if any(strcmp(events,'SaccadeRate'))
            sr = zscore(results.saccades.saccade_rate) + mean(results.diam_left);
            plot(minutes(results.saccades.t(idx_keep) / 60000), ...
                sr(idx_keep), ...
               'Color', saccade_rate_clr, ...
               'LineWidth', 1.5);
        end
    end
    
    % Plot round lines
    if params.plot_rounds && any(strcmp(events,'Rounds'))
        
        % First cycle
        hh = text(minutes(0.05), double(maxx-0.01*maxx),'1.1');
            set(hh,'Color', round_line_color);
            set(hh,'FontWeight', 'bold');
        
        % Cycles
        for j = 1 : length(results.sim2track.cycle_times)
            x1 = results.sim2track.cycle_times(j);
            x1 = minutes(x1 / 60000);
            
            if x1 < t_end
           
                line([x1 x1], [-500 500], ...
                     'Color', round_line_color, ...
                     'LineWidth', 1.5);

                hh = text((x1+0.05), double(maxx-0.01*maxx),sprintf('%d.1',j+1));
                set(hh,'Color', round_line_color);
                set(hh,'FontWeight', 'bold');
             
            end
        end
        
        % Repeats
        cycle = 1; repeat = 2;
        for j = 1 : length(results.sim2track.repeat_times)
           x1 = results.sim2track.repeat_times(j);
           
           while cycle <= length(results.sim2track.cycle_times) && ...
                 results.sim2track.cycle_times(cycle) < x1
               cycle = cycle + 1;
               repeat = 2;
           end
           
           x1 = minutes(x1 / 60000);
           
           if x1 < t_end
               line([x1 x1], [-500 500], ...
                     'Color',round_line_color, ...
                     'LineStyle', '--', ...
                     'LineWidth', 1.5);

               hh = text((x1+0.05), double(maxx-0.01*maxx),sprintf('%d.%d',cycle,repeat));
               set(hh,'Color', round_line_color);
               set(hh,'FontWeight', 'bold');
           end
           
           repeat = repeat + 1;
        end
        
        
    end
    
    % Plot EEG data
    if nargin > 3
        
        N_ch = length(channels);
        clrs = zeros(N_ch, 3);
        
        for d = 1 : length(eeg_data_in)
            eeg_data = eeg_data_in{d};
            
            idx_channels = [];
            for i = 1 : length(channels)
                idx_channels = [idx_channels find(strcmp(eeg_data_in{d}.label,channels{i}))];
            end
            
            N_ch = length(idx_channels);
            
            if params.eeg.interpolate
                L = 0;
                lens = zeros(length(eeg_data.time),1);
                for c = 1 : length(eeg_data.time)
                    lens(c) = length(eeg_data.time{c});
                    L = L + lens(c);
                end
                x_eeg = zeros(L,1);
                y_eeg = zeros(N_ch,L);
                idx = 1;
                for c = 1 : length(eeg_data.time)
                    x_eeg(idx:idx+lens(c)-1) = eeg_data.time{c} / 60;
                    yy = eeg_data.trial{c};
                    if ~isempty(idx_channels)
                        try
                            y_eeg(:,idx:idx+lens(c)-1) = yy(idx_channels,:); 
                        catch wtf
                            a = 0;
                        end
                    end
                    idx = idx+lens(c);
                end
                if params.eeg.smooth > 0
                    for k = 1 : size(y_eeg,1)
                        y_eeg(k,:) = smooth(y_eeg(k,:), params.eeg.smooth, 'moving');
                    end
                end
            else
                L = 0;
                lens = zeros(length(eeg_data.time),1);
                for c = 1 : length(eeg_data.time)
                    lens(c) = length(eeg_data.time{c});
                    L = L + lens(c) + 1;
                end
                x_eeg = zeros(L,1);
                y_eeg = zeros(N_ch,L);
                idx = 1;
                for c = 1 : length(eeg_data.time)
                    x_eeg(idx:idx+lens(c)-1) = eeg_data.time{c} / 60;
                    x_eeg(idx+lens(c)) = x_eeg(idx+lens(c)-1) + 0.00001;
                    yy = eeg_data.trial{c};
                    if ~isempty(idx_channels)
                       y_eeg(:,idx:idx+lens(c)-1) = yy(idx_channels,:); 
                    end
                    y_eeg(:,idx+lens(c)) = nan;
                    idx = idx+lens(c)+1;
                end
            end

            x_eeg = minutes(x_eeg);
            idx_keep2 = x_eeg <= t_end;
            x_eeg = x_eeg(idx_keep2);
            y_eeg = y_eeg(:,idx_keep2);

            itr = params.eeg.height / N_ch;
            gain = itr / params.eeg.stdev; % Std devs per lane

            try
            for i = 1 : N_ch
               y_eeg(i,:) = zscore_nan(y_eeg(i,:)) * gain - i*itr;
            end
            catch err
                fprintf('Error plotting..??');
            end

            hh = plot(x_eeg, y_eeg, 'LineWidth', params.eeg.line_widths(d));
                       
            if d == 1
                ypos = zeros(N_ch,1);
                ystr = cell(N_ch,1);
                for i = 1 : N_ch
                   ypos(i) = nanmean(y_eeg(i,:));
                   ystr(i) = eeg_data.label(idx_channels(i));
                end

                minx = min(y_eeg(:));
                minx = minx-abs(minx)*0.1;

                if ~params.plot_diameter
                    maxx = max(y_eeg(:));
                    maxx = maxx+abs(maxx)*0.1;
                end
                                
                for clr = 1 : N_ch
                    clrs(clr,:) = hh(clr).Color;
                    hh(clr).Color(4) = params.eeg.alpha(d);
                end
                
            else
                for clr = 1 : N_ch
                    hh(clr).Color = clrs(clr,:);
                    hh(clr).Color(4) = params.eeg.alpha(d);
                end
            end
        end
    end
    
    
    % Plot overtake events
    
    if params.plot_overtakes && any(strcmp(events,'Overtakes'))
        
        for j = 1 : length(results.sim2track.overtake_times)
            x1 = results.sim2track.overtake_times(j);
            x1 = minutes(x1 / 60000);
            if x1 < t_end
                line([x1 x1], [-500 500], ...
                     'Color',overtake_color);
            
            end
             
        end
        
    end
    
    
    % Plot saccades
    if params.plot_saccades && any(strcmp(events,'Saccades'))
        for i = 1 : size(results.saccades.saccades,1)
           i_m = results.saccades.saccades(i,2);
           x1 = minutes(results.saccades.t(i_m));
           if x1 < t_end
               line([x1 x1], [-500 500], ...
                     'Color', saccade_line_color, ...
                     'LineWidth', 0.5);
           end
        end
    end

%     grid on;
    
    if params.show_labels
%         hh = title('Pupil diameter with simulation events');
%         set (hh, 'FontSize', 16);
        
        hh = title(plot_title);
        hh.FontSize = 16;

        hh = xlabel('Time (min)');
        set (hh, 'FontSize', 14);

        hh = ylabel('Pupil diameter (mm)');
        set (hh, 'FontSize', 14);
    end
    
    set(gca,'FontName','Arial Narrow');
    set(gca,'FontSize',25);
    set(gca,'tickdir','out');
    set(gca,'ticklength',[.008 .1]);

    
    if params.plot_diameter
        ylim([minx,maxx]);
    end
    
    if nargin > 3
        ylim([minx,maxx]);
        xlim(minutes([0 4]));
        ax = gca;
        
        xlabel('Time (mm:ss)');
        xtickformat('mm:ss');
        ax.XAxis.Label.FontSize = 20;

        ax.YTick=flip(ypos);
        ax.YTickLabel=flip(ystr);
%         ax.XAxis.TickLabelFormat = 'mm:ss'; 
        
    end
    
    if ~isempty(params.xlim)
       xlim(minutes(params.xlim));
    end
    
    resize_window(h, [1500 500], [400 400]);
    
    if ~isempty(params.to_file)
        if iscell(params.to_file)
            for i = 1 : length(params.to_file)
                saveas(h, params.to_file{i});
            end
        else
            saveas(h, params.to_file);
        end
        close(h);
    end

    
    function Z = zscore_nan(X)
        
        Z = (X-nanmean(X)) ./ nanstd(X);
        
    end
    
end

