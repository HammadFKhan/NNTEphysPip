%% Dual Shank LFP
% Plot PA
figure,
subplot(131),imagesc(squeeze(LFP.probe1.PA.Hit)),colormap(myMap),colorbar
subplot(132),imagesc(squeeze(LFP.probe2.PA.Hit)),colormap(myMap),colorbar
subplot(133),imagesc(squeeze(LFP.probe3.PA.Hit)),colormap(myMap),colorbar
%% Plot PA amplitude
time = -1500:1500;
kernal = 75;
supdat1 = smoothdata(mean(squeeze(LFP.probe1.PA.Hit(1:16,:,:))),'gaussian',kernal);
supdat2 = smoothdata(mean(squeeze(LFP.probe3.PA.Hit(1:32,:,:))),'gaussian',kernal);
deepdat1 = smoothdata(mean(squeeze(LFP.probe1.PA.Hit(17:32,:,:))),'gaussian',kernal);
deepdat2 = smoothdata(mean(squeeze(LFP.probe3.PA.Hit(33:64,:,:))),'gaussian',kernal);

figure,
subplot(121),plot(time,supdat1,'k','linewidth', 2),hold on
plot(time,supdat2,'color',[0.7 0.7 0.7],'linewidth', 2)
legend('M2','M1')
box off,set(gca,'tickdir','out','fontsize',12)
xlim([-500 1500]),xlabel('Time from movement (s)'),ylabel('Phase Alignment'),axis square
ylim([0 .5])
subplot(122),plot(time,deepdat1,'k','linewidth', 2),hold on
plot(time,deepdat2,'color',[0.7 0.7 0.7],'linewidth', 2)
legend('M2','M1')
box off,set(gca,'tickdir','out','fontsize',12)
xlim([-500 1500]),xlabel('Time from movement (s)'),ylabel('Phase Alignment'),axis square
ylim([0 .5])
%% Plot PA angle
supdat1 = squeeze(LFP.probe1.PA.Hit_ang(1:16,:,:));
supdat2 = squeeze(LFP.probe3.PA.Hit_ang(1:32,:,:));
deepdat1 = squeeze(LFP.probe1.PA.Hit_ang(17:32,:,:));
deepdat2 = squeeze(LFP.probe3.PA.Hit_ang(33:64,:,:));

figure,
subplot(121),histogram(supdat1,-pi:pi/8:pi,'normalization','probability','edgecolor','none','facecolor',[0 0 0]),hold on
histogram(supdat2,-pi:pi/8:pi,'normalization','probability','edgecolor','none','facecolor',[.7 .7 .7]),set(gca,'TickDir','out'),set(gca,'fontsize',12),box off
xlabel('Max PPA Angle'),ylabel('Probability'),axis square
subplot(122),histogram(deepdat1,-pi:pi/8:pi,'normalization','probability','edgecolor','none','facecolor',[0 0 0]),hold on
histogram(deepdat2,-pi:pi/8:pi,'normalization','probability','edgecolor','none','facecolor',[.7 .7 .7]),set(gca,'TickDir','out'),set(gca,'fontsize',12),box off
xlabel('Max PPA Angle'),ylabel('Probability'),axis square
