function [ h ] = plot_saccades2 ( results, data, params, out2file )

    if nargin < 4
        out2file = false;
    end

    results = results.saccades;
   
    time_sec = results.t / 60000;
    
    gap_clr = [0.95 0.9 0.9];
    saccade_line_color = [0.3 0.3 0.3];
    
    if out2file
        h = figure('visible','off');
    else
        h = figure;
    end
    
    set(h, 'Color', 'w');

    % Plot gaps
    for j = 1 : size(results.tgap,1)
       x1 = time_sec(results.tgap(j,1));
       x2 = x1 + results.tgap(j,2) / 60000;
       if x2 > x1
       rectangle('Position',[x1 -500 x2-x1 1000], ...
                 'EdgeColor','w', ...
                 'FaceColor', gap_clr); 
       end
    end
    
    hold on;
    
    if params.saccades.plot_lines
        for i = 1 : size(results.saccades,1)
           i_m = results.saccades(i,3);
           x1 = time_sec(i_m);
           line([x1 x1], [-500 500], ...
                 'Color', saccade_line_color, ...
                 'LineWidth', 0.5);
            
        end
    end
    
    if params.saccades.plot_saccade_rate
        plot(time_sec,zscore(results.saccade_rate),'m')
    end
    
    % Plot eye position x
    plot(time_sec,zscore(results.x_fixed)+5,'b')

    % Plot eye velocity x
%     dx = results.dx;
%     plot(time_sec,0.3*zscore(dx)-5,'g')
    
    % Plot eye position y
    plot(time_sec,zscore(results.y_fixed)+10,'r')

    % Plot eye velocity
    v = results.velocity;
    plot(time_sec,0.3*zscore(v)-5,'g')
    
    
    
    ylim([-5 15]);
    resize_window(h, [1500 500], [400 400]);
    
    if out2file
        outdir = sprintf('%s/%s/%s/figures', params.root_dir, params.output_dir, data.subject);
        if ~exist(outdir,'dir')
           mkdir(outdir); 
        end
        saveas(h, sprintf('%s/saccades.fig', outdir));
        xlim([0 3]);
        saveas(h, sprintf('%s/saccades.png', outdir));
        close(h);
    end
    

end