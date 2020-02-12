
%%  plot fft

save([dirm, '\results_fft_subjects.mat'])
save([dirm, '\results_fft_stats_subjects.mat'])

figure; hold on
plot([freqsoi(1),freqsoi(end)],[0,0],'k')
p(1) = plot(freqsoi,nanmean(cat(1,tstat_left{:}),1),'b');
p(2)= plot(freqsoi,nanmean(cat(1,tstat_right{:}),1),'r');
p(3)= plot(freqsoi,nanmean(cat(1,tstat_right_lbl{:}),1),'r:');
p(4)= plot(freqsoi,nanmean(cat(1,tstat_rbl_left{:}),1),'b:');
p(5)=  plot(freqsoi,nanmean(cat(1,tstat_pass{:}),1),'g');
xlim([freqsoi(1),freqsoi(end)]);
ylim([-3,3]);
ylabel('Relative power change (T-score)')
xlabel('Frequency (Hz)')
legend(p, {'Left-change', 'Right-change', 'Right rel. to left bl',...
    'Right baseline rel. to left','Passing'}); legend boxoff

figure; hold on;
plot(freqsoi,log10(mean(cat(1,fft_bl_right),1)),':r');
plot(freqsoi,log10(mean(cat(1,fft_right),1)),'r');
plot(freqsoi,log10(mean(cat(1,fft_bl_left),1)),'b:');
plot(freqsoi,log10(mean(cat(1,fft_left),1)),'b');
plot(freqsoi,log10(mean(cat(1,fft_wholepass),1)),'g');
plot(freqsoi,log10(mean(cat(1,fft_wholebl),1)),'g:');
ylabel('Power (10log)')
xlabel('Frequency (Hz)')
legend({'Right baseline','Right', 'Left baseline', 'Left',...
    'Whole pass event','Whole baseline'}); legend boxoff