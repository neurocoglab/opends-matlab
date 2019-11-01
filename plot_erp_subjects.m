function [h, idx_bad] = plot_erp_subjects( summary, channel, tlims, vlims, min_var )

    [h,M_all] = plot_tlocked( summary.erp.left_change.tlocked );
    suptitle('ERPs by subject (all trials)');
    h.Position = [100 700 1500 700];
    
    [h,M_easy] = plot_tlocked( summary.erp.left_change.easy.tlocked );
    suptitle('ERPs by subject (easy trials)');
    h.Position = [80 720 1500 700];
    
    [h,M_diff] = plot_tlocked( summary.erp.left_change.difficult.tlocked );
    suptitle('ERPs by subject (difficult trials)');
    h.Position = [60 740 1500 700];
    
    [h,M_neg] = plot_tlocked( summary.erp.left_change.negative.tlocked );
    suptitle('ERPs by subject (negative outcome trials)');
    h.Position = [40 760 1500 700];
    
    [h,M_pos] = plot_tlocked( summary.erp.left_change.positive.tlocked );
    suptitle('ERPs by subject (positive outcome trials)');
    h.Position = [20 780 1500 700];
    
    tt = summary.erp.left_change.tlocked{1}.time*1000;

    h = plot_difference(tt, M_easy, M_diff, [{'Easy'},{'Difficult'}], 'Overtake onset');
    h.Position = [20 780 1000 500];
    
    h = plot_difference(tt, M_pos, M_neg, [{'Positive'},{'Negative'}], 'Overtake onset');
    h.Position = [20 780 1000 500];
    
    function h = plot_difference( tt, M1, M2, vars, tstr)
        
        h = figure;
        h.Color = 'w';
        
        hh = plot(tt, M1);
        hh.Color = 'r';
        hh.LineWidth = 2;
        hold on;
        
        hh = plot(tt,M2);
        hh.Color = 'b';
        hh.LineWidth = 2;
        
        hh = plot(tt,M2-M1);
        hh.Color = 'k';
        hh.LineWidth = 2;
        
        xlim(tlims*1000);
        xlabel('Time (ms)');
        
        legend(vars);
        
        title(sprintf('%s - %s v. %s [%s]', tstr, vars{1}, vars{2}, channel));
        
    end
    
    
    function [h, M] = plot_tlocked( tlocked )
        h = figure;
        h.Color = 'w';

        
        subplot(3,1,1);
        
        for i = 1 : length(summary.erp.subjects)
            tlock_i = tlocked{i};
            if ~isempty(tlock_i)
                c = find(strcmp(tlock_i.label, channel),1);
                hh = plot(tlock_i.time*1000, tlock_i.avg(c,:));
                hold on;
%                 if mean(tlock_i.avg(c,find(abs(tlock_i.time)<0.01))) < 0
                if max(tlock_i.var(c,:)) < 500
                    hh.Color = 'b';
                    alpha(0.5);
                else
                    hh.Color = 'r';
                    alpha(0.5);
                end
            end

        end

        xlim(tlims*1000);
        title(sprintf('Mean ERP [%s]',channel));
        
        subplot(3,1,2);
        
        idx_v = [find(tlock_i.time>vlims(1),1,'first'):find(tlock_i.time>vlims(2),1,'first')];
        
        M = zeros(length(tlocked{1}.time),1);
        denom = 0;
        
        for i = 1 : length(summary.erp.subjects)
            tlock_i = tlocked{i};
            if ~isempty(tlock_i)
                c = find(strcmp(tlock_i.label, channel),1);
                hh = plot(tlock_i.time*1000, tlock_i.var(c,:));
                hold on;
%                 if mean(tlock_i.avg(c,find(abs(tlock_i.time)<0.01))) < 0
                if max(tlock_i.var(c,idx_v)) < min_var
                    hh.Color = 'b';
                    alpha(0.5);
                    M = M + squeeze(tlock_i.avg(c,:))';
                    denom = denom + 1;
                else
                    fprintf('%s\n',summary.erp.subjects{i});
                    hh.Color = 'r';
                end
            end

        end
        
        xlim(tlims*1000);
        xlabel('Time (ms)');
        title(sprintf('Variance of ERP [%s]',channel));
        
        
        subplot(3,1,3);
        M = M / denom;
        hh = plot(tlock_i.time*1000, M');
        hh.Color = 'k';
        hh.LineWidth = 2;

        xlim(tlims*1000);
        xlabel('Time (ms)');
        title(sprintf('Grand mean ERP [%s]',channel));
        
    end

   

end