
figure,plot(t, X);

F=fft(X);

Fs = 500;            % Sampling frequency                    
T = 1/Fs;     

L = length(X);  

P2 = abs(F/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

figure,plot(f,smooth(P1,1000));



%%

up = envelope(X,100);
% up2 = up*0;
up3 = up2;

for i = 1 : size(up,1)
%     up2(i,:) = smooth(up(i,:),1000);
    [~,pp] = findpeaks(up(i,:));
    idx = [1 pp length(t)];
    up3(i,:) = linterp(t(idx), up(i,idx), t);
end

figure,plot(t, X(1,:));

hold on

% hh = plot(t, up(1,:));
% hh.LineWidth = 2;
% 
% hh = plot(t, up2(1,:));
% hh.LineWidth = 2;
% 
% hh = plot(t, up3(1,:));
% hh.LineWidth = 2;

hh = plot(t, smooth(up3(1,:),200,'loess'));
hh.LineWidth = 2;
