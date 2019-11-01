function [ h ] = plot_saccades ( results, data, params )

    time_sec = results.t / 60000;
    
    gap_clr = [0.95 0.9 0.9];
    
    h = figure;
    set(h, 'Color', 'w');

    % Plot gaps
    for j = 1 : size(results.tgap,1)
       x1 = time_sec(results.tgap(j,1));
       x2 = x1 + results.tgap(j,2) / 60000;
       rectangle('Position',[x1 -500 x2-x1 1000], ...
                 'EdgeColor','w', ...
                 'FaceColor', gap_clr); 
    end
    
    hold on;
    
    % Plot eye position x
    plot(time_sec,zscore(results.x_fixed),'b')

    % Plot eye velocity x
    dx = results.dx;
    plot(time_sec,0.3*zscore(dx)-5,'g')
    
    % Plot eye position y
    plot(time_sec,zscore(results.y_fixed)+10,'r')

    % Plot eye velocity y
    dy = results.dy;
    plot(time_sec,0.3*zscore(dy)+5,'g')
    
    ylim([-10 15]);
    resize_window(h, [1500 500], [400 400]);

end