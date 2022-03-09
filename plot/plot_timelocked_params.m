function [ h ] = plot_timelocked_params( params, stats, labels)

    % Boxplots grouped by parameter type

    h = figure;
    h.Position(3:4) = [800 600];
    
    slopes = [stats.tlock_params{1}.slopes;stats.tlock_params{2}.slopes];
    amplitudes = [stats.tlock_params{1}.amplitudes;stats.tlock_params{2}.amplitudes];
    
    boxplotGroup({slopes,amplitudes},'primaryLabels',labels);
    


end

