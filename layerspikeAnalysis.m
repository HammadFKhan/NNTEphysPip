function output = layerspikeAnalysis(Spikes,LFP)
% Spike and TF analysis of lever data
analyzeLFP = 1;
if isempty(LFP), analyzeLFP = 0;
end

%%% Layer specific spike analysis and statistics

% Assign new variables based on good spikes
hitPSTH = Spikes.PSTH.hit.spks(Spikes.goodSpkComponents);
hitFR = Spikes.PSTH.hit.spkRates(Spikes.goodSpkComponents,:);
missPSTH = Spikes.PSTH.miss.spks(Spikes.goodSpkComponents);
missFR = Spikes.PSTH.miss.spkRates(Spikes.goodSpkComponents,:);
showplot = 0;
if showplot
    y =  [max(hitFR,[],2),max(missFR,[],2)];
    figure,hold on
    x1 = rand(size(y,1),1)/5+1;
    x2 = rand(size(y,1),1)/5+2;
    scatter(x1,y(:,1),'k','filled')
    scatter(x2,y(:,2),'k','filled')
    xlim([0 3])
    plot([x1(:)';x2(:)'], [y(:,1)';y(:,2)'], 'k-')
end

Spk = cell2mat(arrayfun(@(x) x.Clusters(Spikes.goodSpkComponents),Spikes,'UniformOutput',false));
l23Idx = cell2mat(arrayfun(@(x) x.spikeDepth<350,Spk,'UniformOutput',false));
l5Idx = cell2mat(arrayfun(@(x) x.spikeDepth>350,Spk,'UniformOutput',false));

hitFRL23 = hitFR(l23Idx,:);
missFRL23 = missFR(l23Idx,:);
hitFRL5 = hitFR(l5Idx,:);
missFRL5 = missFR(l5Idx,:);

hitPSTHL23 = hitPSTH(l23Idx);
missPSTHL23 = missPSTH(l23Idx);
hitPSTHL5 = hitPSTH(l5Idx);
missPSTHL5 = missPSTH(l5Idx);

data1 = squeeze(LFP.probe1.GED.comp2.missLFP(1,:,:));
params.Fs = 1000; %Fs of LFP and spikes since they have been binned
params.tapers = [5 9];
params.fpass = [5 80];
movingwin = [0.5 0.05];
params.pad = 0;
params.trialave = 1;
niter = 10;
%%% Layer 5 analysis
for n = 1:length(missPSTHL5)
    disp(['Neuron: ' num2str(n)]);
    data2 = missPSTHL5{n}';
    [preCl5(n,:),~,~,~,~,~]=coherencypb(data1(1:1500,:),data2(1:1500,:),params);
    [postCl5(n,:),~,~,~,~,f]=coherencypb(data1(1501:3001,:),data2(1501:3001,:),params);
    temp = data2';
    for nnn = 1:niter
        shufdata2 = [];
        for nn = 1:size(data2,2)
            shufdata2 = vertcat(shufdata2,circshift(temp(nn,:),ceil(3000*rand(1))));
        end
        shufdata2 = shufdata2';
        [shufpretemp(nnn,:),~,~,~,~,~]=coherencypb(data1(1:1500,:),shufdata2(1:1500,:),params);
        [shufposttemp(nnn,:),~,~,~,~,f]=coherencypb(data1(1501:3001,:),shufdata2(1501:3001,:),params);
    end
    shufpreCl5(n,:) = mean(shufpretemp);
    shufpostCl5(n,:) = mean(shufposttemp);
end
zpreCl5 = (preCl5-mean(shufpreCl5,2))./std(shufpreCl5,[],2);
zpostCl5 =(postCl5-mean(shufpostCl5,2))./std(shufpostCl5,[],2);
for n = 1:length(missPSTHL23)
    disp(['Neuron: ' num2str(n)]);
    data2 = missPSTHL23{n}';
    [preCl23(n,:),~,~,~,~,~]=coherencypb(data1(1:1500,:),data2(1:1500,:),params);
    [postCl23(n,:),~,~,~,~,f]=coherencypb(data1(1501:3001,:),data2(1501:3001,:),params);
    temp = data2';
    for nnn = 1:niter
        shufdata2 = [];
        for nn = 1:size(data2,2)
            shufdata2 = vertcat(shufdata2,circshift(temp(nn,:),ceil(3000*rand(1))));
        end
        shufdata2 = shufdata2';
        [shufpretemp(nnn,:),~,~,~,~,~]=coherencypb(data1(1:1500,:),shufdata2(1:1500,:),params);
        [shufposttemp(nnn,:),~,~,~,~,f]=coherencypb(data1(1501:3001,:),shufdata2(1501:3001,:),params);
    end
    shufpreCl23(n,:) = mean(shufpretemp);
    shufpostCl23(n,:) = mean(shufposttemp);
end

zpreCl23 = (preCl23-mean(shufpreCl23,2))./std(shufpreCl23,[],2);
zpostCl23 =(postCl23-mean(shufpostCl23,2))./std(shufpostCl23,[],2);
%%%
if showplot
    figure,
    subplot(221),plot(f,smoothdata(mean(zpreCl23)-3),'k','LineWidth',2),box off, set(gca,'tickdir','out'),set(gca,'fontsize',16);ylim([-2 13])
    subplot(222),plot(f,smoothdata(mean(zpostCl23)-3),'k','LineWidth',2),box off, set(gca,'tickdir','out'),set(gca,'fontsize',16);ylim([-2 13])
    subplot(223),plot(f,smoothdata(mean(zpreCl5)-3),'k','LineWidth',2),box off, set(gca,'tickdir','out'),set(gca,'fontsize',16);ylim([-2 13])
    subplot(224),plot(f,smoothdata(mean(zpostCl5)-3),'k','LineWidth',2),box off, set(gca,'tickdir','out'),set(gca,'fontsize',16);ylim([-2 13])
end
end

function [C,f] = layerspecificCoherence(singleLFP,spkPSTH)
params.tapers = [5 9];
params.fpass = [5 80];
movingwin = [0.5 0.05];
params.pad = 0;
params.trialave = 1;
for n = 1:length(spkPSTH)
    spk = spkPSTH{n}';
    [C(n,:),~,~,~,~,f] = coherencypb(singleLFP,spk,params);
end
end