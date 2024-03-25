
figure,errorbar(1:36,nanmean(hitRateTraining),nanstd(hitRateTraining)/sqrt(3),'.-','Linewidth',2),hold on
errorbar(1:36,nanmean(missRateTraining),nanstd(missRateTraining)/sqrt(3),'.-','Linewidth',2)
errorbar(1:36,nanmean(FARateTraining),nanstd(FARateTraining)/sqrt(3),'.-','Linewidth',2)
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
legend('hitRate','missRate','FARate')
ylabel('Task Rate (outcome/min)')
xlabel('Time (min)')
%%
data = hitRateM2Activate;
figure,errorbar(1:length(data),nanmean(data),nanstd(data)/sqrt(4),'.-','Linewidth',2),hold on
data = missRateM2Activate;
errorbar(1:length(data),nanmean(data),nanstd(data)/sqrt(4),'.-','Linewidth',2)
data = FARateM2Activate;
errorbar(1:length(data),nanmean(data),nanstd(data)/sqrt(4),'.-','Linewidth',2)
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
legend('hitRate','missRate','FARate')
ylabel('Task Rate (outcome/min)')
xlabel('Time (min)')