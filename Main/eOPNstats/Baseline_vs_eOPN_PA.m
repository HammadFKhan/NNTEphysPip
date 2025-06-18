% Statistics of optogenetics across angle/amplitude of PA
%% 
load UCLA_chanmap_fixed.mat
linearProbe = find(s.sorted_probe_wiring(:,2)==20);
PA = Baseline_PA;
%%
figure,
shuf = PA.hitPA_angle(:);
t = squeeze(PA.hitPA_angle);
subplot(331),histogram(t,-pi:pi/8:pi,'normalization','probability'),set(gca,'TickDir','out'),set(gca,'fontsize',12),box off
t = squeeze(PA.missPA_angle);
subplot(334),histogram(t,-pi:pi/8:pi,'normalization','probability'),set(gca,'TickDir','out'),set(gca,'fontsize',12),box off
xlabel('Max PPA Angle'),ylabel('Probability')
[~,r] = max(squeeze(PA.MIFAPA),[],2);
t = squeeze(PA.MIFAPA_angle);
subplot(337),histogram(t,-pi:pi/8:pi,'normalization','probability'),set(gca,'TickDir','out'),set(gca,'fontsize',12),box off
xlabel('Max PPA Angle'),ylabel('Probability')
sz = size(squeeze(PA.hitPA));
dat = mean(squeeze(PA.hitPA(linearProbe,:,1:3001)));
subplot(3,3,[2 3]),plot(1:sz(2),dat);set(gca,'TickDir','out'),set(gca,'fontsize',12),box off
dat = mean(squeeze(PA.missPA(linearProbe,:,1:3001)));
subplot(3,3,[5 6]),plot(1:sz(2),dat);set(gca,'TickDir','out'),set(gca,'fontsize',12),box off
dat = mean(squeeze(PA.MIFAPA(linearProbe,:,1:3001)));
subplot(3,3,[8 9]),plot(1:sz(2),dat);set(gca,'TickDir','out'),set(gca,'fontsize',12),box off


%%
figure,
t = squeeze(eOPN_PA.hitPA_angle);
dat = mean(smoothdata(t(linearProbe,:),'gaussian',5),2);
err = std(t(linearProbe,:),[],2)/sqrt(size(t,2));
subplot(3,2,[1 3 5]),errorbar(dat,1:length(linearProbe),err,'.-','horizontal'),axis tight,set(gca,'TickDir','out'),set(gca,'fontsize',12),box off
set(gca, 'YDir','reverse')
subplot(322),histogram(t(linearProbe(4),:),-pi:pi/8:pi,'normalization','probability'),set(gca,'TickDir','out'),set(gca,'fontsize',12),box off,title('Superficial')
subplot(324),histogram(t(linearProbe(10),:),-pi:pi/8:pi,'normalization','probability'),set(gca,'TickDir','out'),set(gca,'fontsize',12),box off,title('Middle')
subplot(326),histogram(t(linearProbe(20),:),-pi:pi/8:pi,'normalization','probability'),set(gca,'TickDir','out'),set(gca,'fontsize',12),box off,title('Deep')
%% Superficial and Deep PPA response angle for Hits
supdat1 = Baseline_PA.hitPA_angle(linearProbe(1:10),:);
supdat2 = eOPN_PA.hitPA_angle(linearProbe(1:10),:);
deepdat1 = Baseline_PA.hitPA_angle(linearProbe(11:22),:);
deepdat2 = eOPN_PA.hitPA_angle(linearProbe(11:22),:);
figure
subplot(121),histogram(supdat1,-pi:pi/8:pi,'normalization','probability','edgecolor','none'),hold on
histogram(supdat2,-pi:pi/8:pi,'normalization','probability','edgecolor','none'),set(gca,'TickDir','out'),set(gca,'fontsize',12),box off
xlabel('Max PPA Angle'),ylabel('Probability')
[pval, k, K] = circ_kuipertest(supdat1(:), supdat2(:), 60, 0);title(['Superficial ' num2str(pval)]);
subplot(122),histogram(deepdat1,-pi:pi/8:pi,'normalization','probability','edgecolor','none'),hold on
histogram(deepdat2,-pi:pi/8:pi,'normalization','probability','edgecolor','none'),set(gca,'TickDir','out'),set(gca,'fontsize',12),box off
xlabel('Max PPA Angle'),ylabel('Probability')
[pval, k, K] = circ_kuipertest(deepdat1(:), deepdat2(:), 60, 0);title(['Deep ' num2str(pval)]);
%% Superficial and Deep PPA Amplitude for Hits
supdat1 = Baseline_PA.hitPA(linearProbe(1:10),:);
supdat2 = eOPN_PA.hitPA(linearProbe(1:10),:);
deepdat1 = Baseline_PA.hitPA(linearProbe(11:22),:);
deepdat2 = eOPN_PA.hitPA(linearProbe(11:22),:);
statSup = ranksum(max(supdat1,[],2),max(supdat2,[],2));
statDeep = ranksum(max(deepdat1,[],2),max(deepdat2,[],2));
time = (-1500:1500)/1000;
figure
subplot(1,2,1),plot(time,smoothdata(mean(supdat1,1),'movmean',50)),hold on,plot(time,smoothdata(mean(supdat2,1),'movmean',50))
title(['Superficial ' num2str(statSup)]),set(gca,'TickDir','out'),set(gca,'fontsize',12),box off,xlim([-.500 1.500]),ylim([0 0.7])
xlabel('Time (s)'),ylabel('PA Magnitude')
subplot(1,2,2),plot(time,smoothdata(mean(deepdat1,1),'movmean',50)),hold on,plot(time,smoothdata(mean(deepdat2,1),'movmean',50))
title(['Deep ' num2str(statDeep)]),set(gca,'TickDir','out'),set(gca,'fontsize',12),box off,xlim([-.500 1.500]),ylim([0 0.7])
xlabel('Time (s)'),ylabel('PA Magnitude')
figure,
subplot(121),customBoxplot([max(supdat1,[],2),max(supdat2,[],2)]),ylim([0.0 0.3]),ylabel('Max PPA')
title('Superficial'),set(gca,'TickDir','out'),set(gca,'fontsize',12),box off
subplot(122),customBoxplot([max(deepdat1,[],2),max(deepdat2,[],2)]),ylim([0.0 0.3]),ylabel('Max PPA')
title('Deep'),set(gca,'TickDir','out'),set(gca,'fontsize',12),box off










