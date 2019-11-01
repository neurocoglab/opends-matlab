N_s = length(summary.timefreq.right_change.all.tlocked);
X = zeros(N_s,43,1751);
V = zeros(N_s,43,1751);
ch = 36;
idx = [541 1];

for i = 1 : N_s
   
    X(i,:,:) = squeeze(summary.timefreq.right_change.all.tlocked{i}.avg(ch,:,:));
    V(i,:,:) = squeeze(summary.timefreq.right_change.all.tlocked{i}.var(ch,:,:));
    
    
    
end

close all;

Y=X(:,idx(2),idx(1));
figure,plot(Y);
title('Mean power');
xlabel('Subjects');
Y=V(:,idx(2),idx(1));
figure,plot(Y);
title('Variance of power');
xlabel('Subjects');

to_keep = Y < 500;

Y=X(:,idx(2),idx(1));
Y = Y(to_keep);
figure,plot(Y);
title('Mean power (culled)');
xlabel('Subjects');
Y=V(:,idx(2),idx(1));
Y = Y(to_keep);
figure,plot(Y);
title('Variance of power (culled)');
xlabel('Subjects');


figure;
a = ceil(sqrt(length(bldata.trial)));
for i = 1 : length(bldata.trial)
    subplot(a,a,i);
    idx_eeg = bldata.trial{i}(:);
    hist(idx_eeg(~isnan(idx_eeg)),40);
    xlim([-200 200]);
end
suptitle('Raw baselines');

figure;
a = ceil(sqrt(length(bldata.trial)));
for i = 1 : length(bldata.trial)
    subplot(a,a,i);
    idx_eeg = bldata.trial{i}(:);
    hist(zscore(idx_eeg(~isnan(idx_eeg))),40);
    xlim([-20 20]);
end
suptitle('Z-scored baselines');


figure;
a = ceil(sqrt(length(mydata.trial)));
for i = 1 : length(mydata.trial)
    subplot(a,a,i);
    idx_eeg = mydata.trial{i}(:);
    hist(idx_eeg(~isnan(idx_eeg)),40);
    xlim([-200 200]);
end
suptitle('Raw passing');

figure;
a = ceil(sqrt(length(mydata.trial)));
for i = 1 : length(mydata.trial)
    subplot(a,a,i);
    idx_eeg = mydata.trial{i}(:);
    hist(zscore(idx_eeg(~isnan(idx_eeg))),40);
    xlim([-20 20]);
end
suptitle('Z-scored passing');

