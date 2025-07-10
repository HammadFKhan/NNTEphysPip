%% eOPN wave analysis
[basePGD,baseRho] = pgdRho(WavesBaseline);
[eOPNPGD,eOPNRho] = pgdRho(WavesOpto);

%%
t = [];
% figure,customBoxplot(t),box off,set(gca,'tickdir','out','fontsize',14)
% Prepare data for violinplot
t = nan(max([size(baseRho,1) size(eOPNRho,1)]),2);
t(1:size(baseRho,1),1) = baseRho(:,2);
t(1:size(eOPNRho,1),2) = eOPNRho(:,2);
colors = [0.6 0.6 0.6;0.85 0.35 0.1];
figure; % Create a new figure
h_violin = violinplot(t, {'Baseline', 'eOPN'},'ShowWhiskers',false,'ShowBox',false,'MarkerSize',10,'ViolinColor' ,colors);
% --- Axis and Figure Customization ---
ax = gca; % Get current axes handle

ax.TickDir = 'out'; % Ticks point outwards
ax.FontSize = 14; % Font size for tick labels
ax.Box = 'off'; % Turn off the box around the plot

% Y-axis label
ylabel('PGD', 'FontSize', 16);

% Y-axis limits (adjust if needed based on your data range)
% You can let MATLAB auto-set or define explicitly like:
% ax.YLim = [0 1.6]; % Example, adjust based on your data max

% Remove X-axis label (the 'Baseline', 'eOPN' are sufficient from violinplot)
ax.XTickLabel = {'Baseline', 'eOPN'}; % Ensure labels are correctly set
ax.XAxis.FontSize = 14; % Font size for x-axis labels
title(['PGD difference ' num2str(ranksum(baseRho(:,2),eOPNRho(:,2)))])
axis square
%%
figure,errorbar(1:2,[nanmean(baseRho(:,2)) nanmean(eOPNRho(:,2))],[nanstd(baseRho(:,2)) nanstd(eOPNRho(:,2))]),xlim([0.9 2.1])
%% Pre/Post wave Dir
Waves = WavesOpto;
for n = 1:length(Waves.wavesHit)
    idx = Waves.wavesHit(n).evaluationPoints<1500;
    dirWavePreOpto{n} = Waves.wavesHit(n).waveDir(idx);
    idx = Waves.wavesHit(n).evaluationPoints<2500 & Waves.wavesHit(n).evaluationPoints>1501;
    dirWavePostOpto{n} = Waves.wavesHit(n).waveDir(idx);
end
Waves = WavesBaseline;
for n = 1:length(Waves.wavesHit)
    idx = Waves.wavesHit(n).evaluationPoints<1500;
    dirWavePreBase{n} = Waves.wavesHit(n).waveDir(idx);
    idx = Waves.wavesHit(n).evaluationPoints<2500 & Waves.wavesHit(n).evaluationPoints>1501;
    dirWavePostBase{n} = Waves.wavesHit(n).waveDir(idx);
end
%%
figure
histogram(horzcat(dirWavePostBase{:}),-pi:pi/8:pi,'normalization','probability','edgecolor','none'),hold on
histogram(horzcat(dirWavePostOpto{:}),-pi:pi/8:pi,'normalization','probability','edgecolor','none')
box off,set(gca,'tickdir','out','fontsize',16),ylabel('Probability'),xlabel('Phase Angle')
[pval, k, K] = circ_kuipertest(horzcat(dirWavePostBase{:}), horzcat(dirWavePostOpto{:}), 60, 0)

%% Wave speed
Waves = WavesOpto;
for n = 1:length(Waves.wavesHit)
    idx = Waves.wavesHit(n).evaluationPoints<1500;
    speedWavePreOpto{n} = Waves.wavesHit(n).speed(idx);
    idx = Waves.wavesHit(n).evaluationPoints<2500 & Waves.wavesHit(n).evaluationPoints>1501;
    speedWavePostOpto{n} = Waves.wavesHit(n).speed(idx);
end
Waves = WavesBaseline;
for n = 1:length(Waves.wavesHit)
    idx = Waves.wavesHit(n).evaluationPoints<1500;
    speedWavePreBase{n} = Waves.wavesHit(n).speed(idx);
    idx = Waves.wavesHit(n).evaluationPoints<2500 & Waves.wavesHit(n).evaluationPoints>1501;
    speedWavePostBase{n} = Waves.wavesHit(n).speed(idx);
end
%%
figure
histogram(horzcat(speedWavePostBase{:}),'normalization','probability','edgecolor','none'),hold on
histogram(horzcat(speedWavePostOpto{:}),'normalization','probability','edgecolor','none')
box off,set(gca,'tickdir','out','fontsize',16),ylabel('Probability'),xlabel('Phase Angle')
%%
dat1 = horzcat(speedWavePostBase{:});
dat2 = horzcat(speedWavePostOpto{:});
datmean  = [nanmean(dat1) nanmean(dat2)];
datstd = [nanstd(dat1) nanstd(dat2)];
figure,errorbar(1:2,datmean,datstd),xlim([0.9 2.1])
box off,set(gca,'tickdir','out','fontsize',16),ylabel('Wave Speed (cm/s)'),title([num2str(ranksum(dat1,dat2))])

figure,subplot(121),errorbar(1:2,[nanmean(horzcat(speedWavePreBase{:})) nanmean(horzcat(speedWavePostBase{:}))],[nanstd(horzcat(speedWavePreBase{:})) nanstd(horzcat(speedWavePostBase{:}))]),xlim([0.9 2.1])
subplot(122),errorbar(1:2,[nanmean(horzcat(speedWavePreOpto{:})) nanmean(horzcat(speedWavePostOpto{:}))],[nanstd(horzcat(speedWavePreOpto{:})) nanstd(horzcat(speedWavePostOpto{:}))]),xlim([0.9 2.1])
%%
t = [];
t = nan(max([size(dat1,2) size(dat2,2)]),2);
t(1:size(dat1,2),1) = dat1;
t(1:size(dat2,2),2) = dat2;
figure,violinplot(t, {'Baseline', 'eOPN'},'ShowWhiskers',false,'ShowBox',false,'MarkerSize',10,'ViolinColor' ,colors),box off,set(gca,'tickdir','out','fontsize',14)
ylim([0 40])
axis square

%%
figure,customErrorplot(t),box off,set(gca,'tickdir','out','fontsize',14),xlim([0.9 2.1])
%% LOCAL FUNCTIONS

function [normPGDchange,normrhoWave] = pgdRho(Waves)
PGD = arrayfun(@(x) x.PGD, Waves.wavesHit, 'UniformOutput', false);
PGD = vertcat(PGD{:});
rhoWavePre = {};
rhoWavePost = {};
for n = 1:length(Waves.wavesHit)
    idx = Waves.wavesHit(n).evaluationPoints<1500;
    rhoWavePre{n} = cellfun(@mean,Waves.wavesHit(n).rho(idx));
    idx = Waves.wavesHit(n).evaluationPoints<2500 & Waves.wavesHit(n).evaluationPoints>1501;
    rhoWavePost{n} = cellfun(@mean,Waves.wavesHit(n).rho(idx));
end
rhoWave = [cellfun(@mean,rhoWavePre)' cellfun(@mean,rhoWavePost)'];
normrhoWave = [rhoWave(:,1)/mean(rhoWave(:,1))-1 rhoWave(:,2)/mean(rhoWave(:,1))-1]*100;
figure,
errorbar(1:2,nanmean(normrhoWave,1),nanstd(normrhoWave)/sqrt(size(normrhoWave,1)),'-k','linewidth',2),hold on,xlim([0.9 2.1])
for n = 1:size(normrhoWave,1)
    line(1:2,[normrhoWave(n,1) normrhoWave(n,2)]);
end
box off,set(gca,'tickdir','out','fontsize',14),title(['normRhoWave ' num2str(ranksum(normrhoWave(:,1),normrhoWave(:,2)))])

PGDpreCue = mean(PGD(:,1:1500),2);
PGDpostCue = mean(PGD(:,1501:2000),2);
normPGDchange = ([PGDpreCue/mean(PGDpreCue)-1 PGDpostCue/mean(PGDpreCue)-1])*100;
figure,
errorbar(1:2,mean(normPGDchange,1),std(normPGDchange)/sqrt(size(normPGDchange,1)),'-k','linewidth',2),hold on,xlim([0.9 2.1])
for n = 1:length(PGDpreCue)
    line(1:2,[normPGDchange(n,1) normPGDchange(n,2)]);
end
box off,set(gca,'tickdir','out','fontsize',14),title(['normPGD ' num2str(ranksum(normPGDchange(:,1),normPGDchange(:,2)))])
ylabel('PGD change (%)')

end