function [ h ] = plot_saccades_eye ( data, params, out2file )

    if nargin < 3
        out2file = false;
    end
   
    time_sec = data.eye.t / 60000;
    
    gap_clr = [0.95 0.9 0.9];
    saccade_line_color = [0.3 0.3 0.3];
    
    if out2file
        h = figure('visible','off');
    else
        h = figure;
    end
    
    set(h, 'Color', 'w');

    % Plot gaps
    for j = 1 : size(data.eye.tgap,1)
       x1 = time_sec(data.eye.tgap(j,1));
       x2 = x1 + data.eye.tgap(j,2) / 60000;
       if x2 > x1
       rectangle('Position',[x1 -500 x2-x1 1000], ...
                 'EdgeColor','w', ...
                 'FaceColor', gap_clr); 
       end
    end
    
    hold on;
    
    if params.eye.saccades.plots.show_lines
        for i = 1 : size(data.eye.saccades.saccades,1)
           i_m = data.eye.saccades.saccades(i,3);
           x1 = time_sec(i_m);
           line([x1 x1], [-500 500], ...
                 'Color', saccade_line_color, ...
                 'LineWidth', 0.5);
            
        end
    end
    
    if params.eye.saccades.plots.show_rate
        plot(time_sec,zscore(data.eye.saccades.saccade_rate),'m')
    end
    
    % Plot eye position x
    plot(time_sec,zscore(data.eye.saccades.x_fixed)+5,'b')

    % Plot eye velocity x
%     dx = data.eye.saccades.dx;
%     plot(time_sec,0.3*zscore(dx)-5,'g')
    
    % Plot eye position y
    plot(time_sec,zscore(data.eye.saccades.y_fixed)+10,'r')

    % Plot eye velocity
    v = data.eye.saccades.velocity;
    plot(time_sec,0.3*zscore(v)-5,'g')
    
    % Plot threshold
    vt = 0.3*params.eye.saccades.velocity_thres/std(v);
    plot([time_sec(1) time_sec(end)],[vt-5 vt-5],'k--')
    
    
    ylim([-5 15]);
    h.Position = [400 400 1500 500];
    
    xlabel('Time (min)');
  
    if out2file
        outdir = sprintf('%s/%s/figures', params.io.output_dir, data.subject);
        if ~exist(outdir,'dir')
           mkdir(outdir); 
        end
        saveas(h, sprintf('%s/saccades_eye.fig', outdir));
        xlim([0 3]);
        saveas(h, sprintf('%s/saccades_eye.png', outdir));
        close(h);
    end
    

end