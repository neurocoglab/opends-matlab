function h = plot_timelocked_stats( params, stats, groups, layout, sig_dims, sig_clr )
    % Plots time-locked statistics


    plotly_clrs = params.general.plots.plotly_colors;

    x = stats.timelock_stats.time;
    y = stats.timelock_stats.stat;

    h = plotlyfig('Visible','off');
    [data_line,data_shade] = get_plotly_ci_data(x, y, [], [], ...
                                                groups, plotly_clrs, ...
                                                0.1);

    for ii = 1 : length(data_line)
        data_line{ii}.mode='lines';
        data_line{ii}.line.width=2;
    end

    h.data = [data_line;data_shade];

    % Get significance rectangle
    sig_int = stats.timelock_stats.mask;
    diffs = diff([0 sig_int 0]);
    idx_start = find(diffs>0);
    idx_end = find(diffs<0)-1;

    if ~(any(idx_start<1) || any(idx_end<idx_start))
        data_sig = cell(length(idx_start),1);
        if iscell(sig_clr)
            sclr = sig_clr{1}; 
        else
            sclr = sig_clr;
        end
        for i = 1 : length(idx_start)
            xi = [stats.timelock_stats.time(idx_start(i)), ...
                  stats.timelock_stats.time(idx_end(i))]';

            yi = [(sig_dims(1)+sig_dims(2)/2);(sig_dims(1)+sig_dims(2)/2); ...
                  (sig_dims(1)-sig_dims(2)/2);(sig_dims(1)-sig_dims(2)/2)];
            X = struct(...
                                    'x', [xi;flip(xi)], ...
                                    'y', yi, ...
                                    'name', ['sig_' i], ...
                                    'visible', 1, ...
                                    'fill', 'toself', ...
                                    'mode', 'lines', ...
                                    'fillcolor', sclr, ...
                                    'line', struct('color', sclr, 'width', 0), ...
                                    'showlegend', false ...
                                    );
            data_sig(i) = {X};
        end
        h.data = [data_sig;h.data];
    end

    if isfield(stats.timelock_stats, 'mask2')
        % Get lesser significance rectangle
        if iscell(sig_clr)
            sclr = sig_clr{2}; 
        else
            sclr = sig_clr;
        end
        sig_int = stats.timelock_stats.mask2;
        diffs = diff([0 sig_int 0]);
        idx_start = find(diffs>0);
        idx_end = find(diffs<0)-1;

        if ~(any(idx_start<1) || any(idx_end<idx_start))
            data_sig = cell(length(idx_start),1);
            for i = 1 : length(idx_start)
                xi = [stats.timelock_stats.time(idx_start(i)), ...
                      stats.timelock_stats.time(idx_end(i))]';

                yi = [(sig_dims(1)+sig_dims(2)/2);(sig_dims(1)+sig_dims(2)/2); ...
                      (sig_dims(1)-sig_dims(2)/2);(sig_dims(1)-sig_dims(2)/2)];
                X = struct(...
                                        'x', [xi;flip(xi)], ...
                                        'y', yi, ...
                                        'name', ['sig_' i], ...
                                        'visible', 1, ...
                                        'fill', 'toself', ...
                                        'mode', 'lines', ...
                                        'fillcolor', sclr, ...
                                        'line', struct('color', sclr, 'width', 0), ...
                                        'showlegend', false ...
                                        );
                data_sig(i) = {X};
            end
            h.data = [data_sig;h.data];
        end
    end
    
    h.layout = layout;
    h.layout.yaxis.title = stats.statistic;
    h.layout.xaxis.title = 'Time relative to button press (s)';
        
        
    end
