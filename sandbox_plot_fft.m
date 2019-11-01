figure

for i = 1 : length(data_epochs.baseline.time)
   
    t = data_epochs.baseline.time{i};
    X = data_epochs.baseline.trial{i};
    X = X(20,:);
    plot(t,X);
    hold on;
    
end