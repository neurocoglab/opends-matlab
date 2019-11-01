close all;

[H,X] = hist(stats.left_change.timelockstats.posdistribution,40);
Y = H / trapz(H);
C=zeros(40,1);
for i = 1 : 40
    C(i) = trapz(Y(1:i));
end

figure, plot(X,Y);
hold on;
plot(X,C);

idx_crit = find(C>=0.975,1);
v_crit = X(idx_crit);

hold on;
line([v_crit v_crit],[0 1],'Color','g');

fprintf('Critical stat (97.5%%) is %1.2f\n', v_crit);