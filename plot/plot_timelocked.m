function h = plot_timelocked( params, stats, groups, layout, sig_dims, sig_clr )
    
    plotly_clrs = params.general.plots.plotly_colors;
    N = length(stats.subjects);
%     plotly_layouts = load(params.general.plots.plotly_layouts_file);

    if length(stats.timelock_avr) == 1
        x = stats.timelock_avr{1}.time;
        y = stats.timelock_avr{1}.avg;
        var = stats.timelock_avr{1}.var;
        lower = y - 1.96*sqrt(var/N);
        upper = y + 1.96*sqrt(var/N);
    else
        N_s = length(stats.timelock_avr);
        N_t = length(stats.timelock_avr{1}.time);
        x = zeros(N_t,N_s);
        y = zeros(N_t,N_s);
        var = zeros(N_t,N_s);
        upper = zeros(N_t,N_s);
        lower = zeros(N_t,N_s);
        for ii = 1 : length(stats.timelock_avr)
            x(:,ii) = stats.timelock_avr{ii}.time;
            y(:,ii) = stats.timelock_avr{ii}.avg;
            var(:,ii) = stats.timelock_avr{ii}.var;
            lower(:,ii) = y(:,ii) - 1.96*sqrt(var(:,ii)/N);
            upper(:,ii) = y(:,ii) + 1.96*sqrt(var(:,ii)/N);
        end
    end

    h = plotlyfig('Visible','off');
    [data_line,data_shade] = get_plotly_ci_data(x, y, lower, upper, ...
                                                groups, plotly_clrs, ...
                                                0.1);

    for ii = 1 : length(data_line)
        data_line{ii}.mode='lines';
        data_line{ii}.line.width=2;
    end

    h.data = [data_line;data_shade];

    % Get significance rectangle
    sig_int = stats.timelockstats.mask;
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
            xi = [stats.timelockstats.time(idx_start(i)), ...
                  stats.timelockstats.time(idx_end(i))]';

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

    if isfield(stats.timelockstats, 'mask2')
        % Get lesser significance rectangle
        if iscell(sig_clr)
            sclr = sig_clr{2}; 
        else
            sclr = sig_clr;
        end
        sig_int = stats.timelockstats.mask2;
        diffs = diff([0 sig_int 0]);
        idx_start = find(diffs>0);
        idx_end = find(diffs<0)-1;

        if ~(any(idx_start<1) || any(idx_end<idx_start))
            data_sig = cell(length(idx_start),1);
            for i = 1 : length(idx_start)
                xi = [stats.timelockstats.time(idx_start(i)), ...
                      stats.timelockstats.time(idx_end(i))]';

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
    h.layout.yaxis.title = 'PD (normalized to baseline)';
    h.layout.xaxis.title = 'Time relative to button press (s)';
        
        
    end
