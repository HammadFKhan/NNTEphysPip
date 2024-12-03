%%
colors = [0 0.4470 0.7410;0.75 0.75 0.75;190/255 30/255 45/255];
time = -1499:20:1500;
figure,plot(time(2:end),M1neuralDynamics.neuralDiffhitmiss(:,1:2),'color',colors(1,:),'LineWidth',2),ylim([-0.02 0.04]), hold on
plot(time(2:end),M1neuralDynamics.neuralDiffMI(:,1:2),'color',colors(2,:),'LineWidth',2),ylim([-0.02 0.04])
box off,set(gca,'tickdir','out','fontsize',14),axis square,xlim([-500 1500])
ylabel('Sequence Distance'),xlabel('Time (s)')

colors = [0 0.4470 0.7410;0.75 0.75 0.75;190/255 30/255 45/255];
time = -1499:20:1500;
figure,plot(time(2:end),M2neuralDynamics.neuralDiffhitmiss(:,1:2),'color',colors(1,:),'LineWidth',2),ylim([-0.06 0.16]), hold on
plot(time(2:end),M2neuralDynamics.neuralDiffMI(:,1:2),'color',colors(2,:),'LineWidth',2),ylim([-0.06 0.16])
box off,set(gca,'tickdir','out','fontsize',14),axis square,xlim([-500 1500])
ylabel('Sequence Distance'),xlabel('Time (s)')
%% Plot neural trajectory speed in relation to RT
neuralDynamics = M2neuralDynamics;
figure,scatter(neuralDynamics.PQ,neuralDynamics.rawreactionTime,'k','filled')
xlim([0 15])
ylim([0 2])
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
xlabel('Neural trajectory speed'),ylabel('Reaction time (s)')
mdl = fitlm(neuralDynamics.PQ,neuralDynamics.rawreactionTime)
figure,plot(mdl)
xlim([0 15])
ylim([0 2])
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
xlabel('Neural trajectory speed'),ylabel('Reaction time (s)')
legend(num2str(mdl.Rsquared.Ordinary),num2str(mdl.Coefficients.pValue(2)))
title('')
%% Plot neural trajectory STD
neuralDynamics = M1neuralDynamics;
temp = nan(max([size(neuralDynamics.hit.s,2),size(neuralDynamics.miss.s,2),size(neuralDynamics.MIFA.s,2)]),3);
temp(1:size(neuralDynamics.hit.s,2),1) = neuralDynamics.hit.s(1,:);
temp(1:size(neuralDynamics.miss.s,2),2) = neuralDynamics.miss.s(1,:);
temp(1:size(neuralDynamics.MIFA.s,2),3) = neuralDynamics.MIFA.s(1,:);

figure,customBarplot(temp);
box off,set(gca,'tickdir','out','fontsize',14),axis square
ylabel('Trajectory Deviation')
[p,t,stats] = anova1(temp)
c = multcompare(stats)