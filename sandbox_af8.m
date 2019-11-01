close all;

channels = [{'AF7'},{'AF8'}];
N_subj = length(subjects);
N_ch = length(channels);
N_t = length(summary.erp.left_change.tlocked{1}.time);

values = zeros(N_subj,N_ch,N_t);
vars = zeros(N_subj,N_ch,N_t);

for i = 1 : N_subj
    
    tl = summary.erp.left_change.tlocked{i};
    
    for c = 1 : N_ch
        idx_c = find(strcmp(tl.label, channels{c}));
        values(i,c,:) = tl.avg(idx_c,:);
        vars(i,c,:) = tl.var(idx_c,:);
    end
    
end

for c = 1 : N_ch
    figure,plot(summary.erp.left_change.tlocked{1}.time,squeeze(values(:,c,:)));
    title(channels{c});
    
    figure,plot(mean(squeeze(values(:,c,:)),2));
    hold on;
    plot(mean(squeeze(vars(:,c,:)),2));
    title(['By subjects: ' channels{c}]);
        
end

badsubs = [6 17];

for s = 1 : length(badsubs)
    figure;
    tl = summary.erp.left_change.tlocked{badsubs(s)};
    for c = 1 : N_ch
        plot(summary.erp.left_change.tlocked{1}.time,squeeze(values(badsubs(s),c,:)));
        hold on;
    end
    
    legend(channels);
    title(['Mean signal: ' subjects{badsubs(s)}]);
    
    figure;
    tl = summary.erp.left_change.tlocked{badsubs(s)};
    for c = 1 : N_ch
        plot(summary.erp.left_change.tlocked{1}.time,squeeze(vars(badsubs(s),c,:)));
        hold on;
    end
    
    legend(channels);
    title(['Signal variance: ' subjects{badsubs(s)}]);
    
end
