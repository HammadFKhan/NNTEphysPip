%% Sequence analysis for supplementary figure
% Taken from: Neural Sequences as an Optimal Dynamical Regime for the Readout of Time
% https://doi.org/10.1016/j.neuron.2020.08.020
% To quantify the degree to which different dynamical regimes could be
% described as a sequence we developed a novel measure of sequentiality
% index (SqI). Our SqI was inspired by a previously described measure
% (Orhan and Ma, 2019), but relies on temporal sparsity rather than the
% ridge-to-background ratio, and is bounded between 0 and 1 (Figures 2 and
% 4). As illustrated in Figure 2D, this normalized sequentiality index was
% calculated as: PE = XM j = 1
%  pj logðpjÞ
% , logðMÞ TS = 1  < XN i = 1
%  r
% t i logðr t iÞ , logðNÞ>t SqI =
% ffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffi PE  TS p where PE
% (peak entropy) captures the entropy of the distribution of peak times
% across the entire population. M is the number of bins for estimating the
% peak time distribution (30 and 60 for neural dynamics in the short and
% long intervals in Figure 2, and 10 for the simulated dynamics in Figure
% 4), and pj is the number of units that peak in time bin j divided by the
% total number of units. TS (temporal sparsity) provides a measure of
% entropy of the distribution of normalized activity in any given bin. TS
% is maximum when only one unit accounts for all the activity in each time
% bin t. rt i represents the activity of unit i in time bin t, normalized
% by the sum of the activity of all neurons at time t;<> t denotes the time
% average; and N is the number of units. The SqI approaches 1 when the peak
% times of each neuron homogeneously tile the entire duration and one
% neuron is active at every moment in time withthe temporal fields being
% non-overlapping. Note that SqI is sensitive to the width of the temporal
% fields of the neurons, and this is an important feature because the
% extreme sequences with very broad temporal field converge toward a regime
% with a fixed-point attractor. We stress that because sequentiality
% represents a fairly complex feature of dynamics, that no single measure
% of sequentiality is perfect. The measure defined here seems effective at
% capturing the biological dynamical regimes studied here. However, it is
% important to note that this measure can break down under some artificial
% manipulations (e.g., if one were to shuffle the time bins of a sequential
% neural trajectory).
% Here we use Bin width of 20
load('D:\SequenceDat\M1Sq.mat')
load('D:\SequenceDat\M2Sq.mat')
colors = [0 0.4470 0.7410;0.75 0.75 0.75;0 0.4470 0.7410;190/255 30/255 45/255]; % data is indexed as Hit, Miss, MIHIT, MIFA

dat = [];
for n = 1:length(M1SqEntropyMice(1).SqEntropy.CueHit.bin)
    bin = n;
datM1 = [arrayfun(@(x) x.SqEntropy.CueHit.SqI(bin), M1SqEntropyMice)',arrayfun(@(x) x.SqEntropy.CueMiss.SqI(bin), M1SqEntropyMice)',arrayfun(@(x) x.SqEntropy.MIHit.SqI(bin), M1SqEntropyMice)', arrayfun(@(x) x.SqEntropy.MIFA.SqI(bin), M1SqEntropyMice)'];
dat(n,:) = mean(datM1,1);
end
figure;
for n = 1:size(dat,2)
    plot(M1SqEntropyMice(1).SqEntropy.CueHit.bin,dat(:,n),'color',colors(n,:));hold on
end

 box off, set(gca, 'tickdir','out','fontsize',10),axis square


dat = [];
for n = 1:length(M2SqEntropyMice(1).SqEntropy.CueHit.bin)
    bin = n;
datM2 = [arrayfun(@(x) x.SqEntropy.CueHit.SqI(bin), M2SqEntropyMice)',arrayfun(@(x) x.SqEntropy.CueMiss.SqI(bin), M2SqEntropyMice)',arrayfun(@(x) x.SqEntropy.MIHit.SqI(bin), M2SqEntropyMice)', arrayfun(@(x) x.SqEntropy.MIFA.SqI(bin), M2SqEntropyMice)'];
dat(n,:) = mean(datM2,1);
end
figure;
for n = 1:size(dat,2)
    plot(M2SqEntropyMice(1).SqEntropy.CueHit.bin,dat(:,n),'color',colors(n,:));hold on
end
 box off, set(gca, 'tickdir','out','fontsize',10),axis square
%% Sequentiality Index
bin = 1;
datM1 = [arrayfun(@(x) x.SqEntropy.CueHit.SqI(bin), M1SqEntropyMice)',arrayfun(@(x) x.SqEntropy.CueMiss.SqI(bin), M1SqEntropyMice)',arrayfun(@(x) x.SqEntropy.MIHit.SqI(bin), M1SqEntropyMice)', arrayfun(@(x) x.SqEntropy.MIFA.SqI(bin), M1SqEntropyMice)'];
bin = 15;
datM2 = [arrayfun(@(x) x.SqEntropy.CueHit.SqI(bin), M2SqEntropyMice)',arrayfun(@(x) x.SqEntropy.CueMiss.SqI(bin), M2SqEntropyMice)',arrayfun(@(x) x.SqEntropy.MIHit.SqI(bin), M2SqEntropyMice)', arrayfun(@(x) x.SqEntropy.MIFA.SqI(bin), M2SqEntropyMice)'];

f = figure;
for n = 1:4
    subplot(1,4,n)
    errorbar(1,mean(datM1(:,n)),std(datM1(:,n)),'ko'),hold on
    scatter(1*ones(size(datM1,1),1),datM1(:,n),10,colors(n,:),'filled','jitter','on','jitterAmount',0.1)
    
    errorbar(2,mean(datM2(:,n)),std(datM2(:,n)),'ko'),hold on
    scatter(2*ones(size(datM2,1),1),datM2(:,n),10,colors(n,:),'filled','jitter','on','jitterAmount',0.1)
    xlim([0.5 2.5])
    box off, set(gca, 'tickdir','out','fontsize',10),axis square
    ylabel('Sequentiality Index')
end
f.Position = [681 559 860 220];

%% Temporal Sparsity
bin = 1;
datM1 = [arrayfun(@(x) x.SqEntropy.CueHit.TS(bin), M1SqEntropyMice)',arrayfun(@(x) x.SqEntropy.CueMiss.TS(bin), M1SqEntropyMice)',arrayfun(@(x) x.SqEntropy.MIHit.TS(bin), M1SqEntropyMice)', arrayfun(@(x) x.SqEntropy.MIFA.TS(bin), M1SqEntropyMice)'];
bin = 15;
datM2 = [arrayfun(@(x) x.SqEntropy.CueHit.TS(bin), M2SqEntropyMice)',arrayfun(@(x) x.SqEntropy.CueMiss.TS(bin), M2SqEntropyMice)',arrayfun(@(x) x.SqEntropy.MIHit.TS(bin), M2SqEntropyMice)', arrayfun(@(x) x.SqEntropy.MIFA.TS(bin), M2SqEntropyMice)'];

f = figure;
for n = 1:4
    subplot(1,4,n)
    errorbar(1,mean(datM1(:,n)),std(datM1(:,n)),'ko'),hold on
    scatter(1*ones(size(datM1,1),1),datM1(:,n),10,colors(n,:),'filled','jitter','on','jitterAmount',0.1)
    
    errorbar(2,mean(datM2(:,n)),std(datM2(:,n)),'ko'),hold on
    scatter(2*ones(size(datM2,1),1),datM2(:,n),10,colors(n,:),'filled','jitter','on','jitterAmount',0.1)
    xlim([0.5 2.5])
    box off, set(gca, 'tickdir','out','fontsize',10),axis square
    ylabel('Temporal Sparsity')
end
f.Position = [681 559 860 220];
%% Peak Entropy
bin = 1;
datM1 = [arrayfun(@(x) x.SqEntropy.CueHit.PE(bin), M1SqEntropyMice)',arrayfun(@(x) x.SqEntropy.CueMiss.PE(bin), M1SqEntropyMice)',arrayfun(@(x) x.SqEntropy.MIHit.PE(bin), M1SqEntropyMice)', arrayfun(@(x) x.SqEntropy.MIFA.PE(bin), M1SqEntropyMice)'];
bin = 15;
datM2 = [arrayfun(@(x) x.SqEntropy.CueHit.PE(bin), M2SqEntropyMice)',arrayfun(@(x) x.SqEntropy.CueMiss.PE(bin), M2SqEntropyMice)',arrayfun(@(x) x.SqEntropy.MIHit.PE(bin), M2SqEntropyMice)', arrayfun(@(x) x.SqEntropy.MIFA.PE(bin), M2SqEntropyMice)'];

f = figure;
for n = 1:4
    subplot(1,4,n)
    errorbar(1,mean(datM1(:,n)),std(datM1(:,n)),'ko'),hold on
    scatter(1*ones(size(datM1,1),1),datM1(:,n),10,colors(n,:),'filled','jitter','on','jitterAmount',0.1)
    
    errorbar(2,mean(datM2(:,n)),std(datM2(:,n)),'ko'),hold on
    scatter(2*ones(size(datM2,1),1),datM2(:,n),10,colors(n,:),'filled','jitter','on','jitterAmount',0.1)
    xlim([0.5 2.5])
    box off, set(gca, 'tickdir','out','fontsize',10),axis square
    ylabel('Peak Entropy')
end
f.Position = [681 559 860 220];
%%
[~,~,stats] = anova1([datM1,datM2]);
c = multcompare(stats)
%% Make gaussian for entropy figure
figure
L=20;
V = linspace(0,L,301);
sigma = L/20;
T0 = exp(-((V-mean(V))/sigma).^2/2);
plot(V,T0);axis off

%% Spike precision by Neurons
load('D:\SequenceDat\M1SpikePrecisionbyNeuron.mat')
load('D:\SequenceDat\M2SpikePrecisionbyNeuron.mat')
colors = [0 0.4470 0.7410;0.75 0.75 0.75;190/255 30/255 45/255]; % data is indexed as Hit, Miss, MIFA
plotOn = 0;
totalSpikes = M1totalSpikes; % Specify brain region here
dat = totalSpikes.hit;
precisionHitSpk = calcSpkPrecision(dat,plotOn);
dat = totalSpikes.miss;
precisionMissSpk = calcSpkPrecision(dat,plotOn);
dat = totalSpikes.FA;
precisionFASpk = calcSpkPrecision(dat,plotOn);
f = figure;
subplot(131),imagesc(precisionHitSpk),caxis([0 1]),axis square, colorbar
subplot(132),imagesc(precisionMissSpk),caxis([0 1]),axis square, colorbar
subplot(133),imagesc(precisionFASpk),caxis([0 1]),axis square, colorbar
f.Position = [681 559 860 220];
figure;
histogram(precisionHitSpk,0:0.01:1,'Normalization','Probability','FaceColor',colors(1,:),'EdgeColor','none'),hold on
histogram(precisionMissSpk,0:0.01:1,'Normalization','Probability','FaceColor',colors(2,:),'EdgeColor','none')
histogram(precisionFASpk,0:0.01:1,'Normalization','Probability','FaceColor',colors(3,:),'EdgeColor','none')
box off, set(gca,'tickdir','out','fontsize',14),axis square
xlabel('Spike Precision'),ylabel('Probability')
%%
figure,customBarplot([nanmean(precisionHitSpk,2),nanmean(precisionMissSpk,2),nanmean(precisionFASpk,2)])
ylabel('Spike Precision')
box off, set(gca,'tickdir','out','fontsize',14)
ylim([0 0.5])
%%% STATS
[~,~,stats] = anova1([nanmean(precisionHitSpk,2),nanmean(precisionMissSpk,2),nanmean(precisionFASpk,2)])
c = multcompare(stats)
%% Within neuron precision
% Since the precision of neurons is done within mice, we are basically
% saying what is the precision of the population response. Ie. neuronal
% sequences that self organize will have higher precision. 
% In this analysis we are actually calculating the precision of a single
% neuron to change as a function of task trial.

totalSpikes = M2totalSpikes; % Specify brain region here
interSpkPrecision = [];
% We combined the spikes but we also need to index so we can find it again
id = 1:size(totalSpikes.hit,1);
combineTaskSpks = [totalSpikes.hit;totalSpikes.miss;totalSpikes.FA]; 
assert(id(end)*3==size(combineTaskSpks,1))
allspkPrecision = calcSpkPrecision(combineTaskSpks,0);
for n = id
    interSpkPrecision(n,:) = [1,allspkPrecision(n,n+id(end)),allspkPrecision(n,n+id(end)*2)];
end
% Fix NaN scenarios because we wrote it in the function as NaN == 1
interSpkPrecision(isnan(interSpkPrecision)) = 1;

dat = interSpkPrecision;
figure,customBarplot(dat,'Scatter','off'),hold on
for n = 1:size(dat,1)
    scatter(1:3,[dat(n,1),dat(n,2),dat(n,3)],'k','filled')
    line([1 2 3],[dat(n,1),dat(n,2),dat(n,3)])
end

ylabel('Spike Precision')
box off, set(gca,'tickdir','out','fontsize',14)
ylim([0 1.1])
%% Plot out example neurons for visualization
colors = [0 0.4470 0.7410;0.75 0.75 0.75;190/255 30/255 45/255];
figure;
neuronId = 120;
time = -1.5:0.001:1.5;
plot(time,mean(combineTaskSpks(1:id(end),:)),'color',colors(1,:),'LineWidth',2);hold on
plot(time,mean(combineTaskSpks(id(end)+1:id(end)*2,:)),'color',colors(2,:),'LineWidth',2);hold on
plot(time,mean(combineTaskSpks(id(end)*2+1:id(end)*3,:)),'color',colors(3,:),'LineWidth',2)
xlim([-0.5 1.5]),axis square,box off,set(gca,'tickdir','out','fontsize',14)
ylim([-0.5 1])



%% LOCAL FUNCTIONS
function precisionSpk = calcSpkPrecision(dat,plotOn)
precisionSpk = abs(corrcoef(dat'));
if plotOn
    figure,imagesc(precisionSpk),caxis([0 1]),axis square, colorbar
    precisionSpk(precisionSpk==1) = NaN;
    figure,histogram(precisionSpk,0:0.01:1,'Normalization','Probability','EdgeColor','none')
    box off, set(gca,'tickdir','out','fontsize',14),axis square
    xlabel('Spike Precision'),ylabel('Probability')
end

end