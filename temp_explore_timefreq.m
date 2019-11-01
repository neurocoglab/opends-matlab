trialdefs = fieldnames(summary.timefreq.left_change);

for i = 1 : length(trialdefs)

    for j = 1 : length(subjects)
        
        tlocked = summary.timefreq.left_change.(trialdefs{i}).tlocked{j};
        
        if isempty(tlocked)
           fprintf('Not enough trials for subject %s, of type "%s"\n', subjects{j}, trialdefs{i});
        end
        
    end

end