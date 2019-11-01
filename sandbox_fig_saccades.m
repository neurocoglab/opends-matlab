h = figure;

oclr = [0.8242 0.6406 0.1719];
saccade_line_color = [0.3 0.3 0.3];

for j = 1 : size(data.eye.tgap,1)
   x1 = results.t(data.eye.tgap(j,1));
   x2 = x1 + data.eye.tgap(j,2);
   hr = rectangle('Position',[x1 -500 x2-x1 1000], ...
                  'EdgeColor','w', ...
                  'FaceColor', gap_clr); 
end

hold on
plot(results.t,zscore(results.saccades.dy),'Color',[0.8 0.8 0.8])
plot(results.t,zscore(results.saccades.x_fixed),'b')
plot(results.t,zscore(results.saccades.y_fixed),'r')
plot(results.t,results.blinks*5,'m')
plot(results.t,-results.blinks*5,'m')

plot(results.t,results.saccades.velocity/40 - 5,'Color',oclr)

for i = 1 : size(results.saccades.saccades,1)
    idx = results.saccades.saccades(i,3);
    x1 = results.t(idx);
           line([x1 x1], [-500 500], ...
                 'Color', saccade_line_color, ...
                 'LineWidth', 0.5);
end

ylim([-5 5]);

resize_window(h,[1300 600],[200 800]);

xlim([0 60000]);