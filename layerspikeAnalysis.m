function Spikes = layerspikeAnalysis(Spikes,Behaviour,LFP)
% Spike and TF analysis of lever data
analyzeLFP = 1;
if isempty(LFP), analyzeLFP = 0;
end

%%% Behaviour parsing
winBeforeCue = 1:Behaviour.parameters.windowBeforeCue*Behaviour.parameters.Fs;
winAfterCue = winBeforeCue(end)+1:Behaviour.parameters.windowAfterCue*Behaviour.parameters.Fs+winBeforeCue(end)+1;
winRT = arrayfun(@(x) winAfterCue(1)+x.reactionTime*Behaviour.parameters.Fs,Behaviour.cueHitTrace);
winRW = arrayfun(@(x) find(x.LFPIndex==x.rewardIndex),Behaviour.cueHitTrace,'UniformOutput',false);
checkRW = cellfun(@(x) isempty(x),winRW);
winRW(checkRW) = cell({0});% Time from cue to Reward
winRW = cell2mat(winRW);
winRW(winRW==0) = NaN;
%%% Layer specific spike analysis and statistics
% Assign new variables based on good spikes
if isfield(Spikes.PSTH,'hit')
hitPSTH = Spikes.PSTH.hit.spks(Spikes.goodSpkComponents);
hitFR = Spikes.PSTH.hit.spkRates(Spikes.goodSpkComponents,:);
end
if isfield(Spikes.PSTH,'miss')
missPSTH = Spikes.PSTH.miss.spks(Spikes.goodSpkComponents);
missFR = Spikes.PSTH.miss.spkRates(Spikes.goodSpkComponents,:);
end
if isfield(Spikes.PSTH,'MIHit')
    mihitPSTH = Spikes.PSTH.MIHit.spks(Spikes.goodSpkComponents);
    mihitFR = Spikes.PSTH.MIHit.spkRates(Spikes.goodSpkComponents,:);
end
if isfield(Spikes.PSTH,'MIFA')
    mifaPSTH = Spikes.PSTH.MIFA.spks(Spikes.goodSpkComponents);
    mifaFR = Spikes.PSTH.MIFA.spkRates(Spikes.goodSpkComponents,:);
end


%%%
Spk = cell2mat(arrayfun(@(x) x.Clusters(Spikes.goodSpkComponents),Spikes,'UniformOutput',false));
l23Idx = cell2mat(arrayfun(@(x) x.spikeDepth<350,Spk,'UniformOutput',false));
l5Idx = cell2mat(arrayfun(@(x) x.spikeDepth>350,Spk,'UniformOutput',false));


hitFRL23 = hitFR(l23Idx,:);
missFRL23 = missFR(l23Idx,:);
mihitFRL23 = mihitFR(l23Idx,:);
mifaFRL23 = mifaFR(l23Idx,:);

hitFRL5 = hitFR(l5Idx,:);
missFRL5 = missFR(l5Idx,:);
mihitFRL5 = mihitFR(l5Idx,:);
mifaFRL5 = mifaFR(l5Idx,:);
% Do for PSTH
hitPSTHL23 = hitPSTH(l23Idx);
missPSTHL23 = missPSTH(l23Idx);
mihitPSTHL23 = mihitPSTH(l23Idx);
mifaPSTHL23 = mifaPSTH(l23Idx);

hitPSTHL5 = hitPSTH(l5Idx);
missPSTHL5 = missPSTH(l5Idx);
mihitPSTHL5 = mihitPSTH(l5Idx);
mifaPSTHL5 = mifaPSTH(l5Idx);

%%% tuning curve characterization (preferred stimulas, supressed, activate
%%% after reward, non-specific. Defined under task specific

hitspkTuning = spkResponse(hitFR,winBeforeCue,winAfterCue,winRT,winRW);
missspkTuning = spkResponse(missFR,winBeforeCue,winAfterCue,winRT,winRW);

%%% output some stuff

Spikes.spikeProp.hitspkTuning = hitspkTuning;
Spikes.spikeProp.missspkTuning = missspkTuning;
Spikes.spikeProp.hitFRL23 = hitFRL23;
Spikes.spikeProp.missFRL23 = missFRL23;
Spikes.spikeProp.hitFRL5 = hitFRL5;
Spikes.spikeProp.missFRL5 = missFRL5;

showplot = 1;
if showplot
    y =  [mean(hitFRL23(:,winBeforeCue(end):winBeforeCue(end)+1000),2),mean(missFRL23(:,winBeforeCue(end):winBeforeCue(end)+1000),2)];
    [~,p] = ttest(y(:,1),y(:,2));
    disp(['L23 Hit vs Miss p = ' num2str(p)])
    figure,subplot(121),hold on
    x1 = rand(size(y,1),1)/5+1;
    x2 = rand(size(y,1),1)/5+2;
    scatter(x1,y(:,1),'r','filled')
    scatter(x2,y(:,2),'r','filled')
    xlim([.75 2.25])
    plot([x1(:)';x2(:)'], [y(:,1)';y(:,2)'], 'k-')
    y =  [mean(hitFRL5(:,winBeforeCue(end):winBeforeCue(end)+1000),2),mean(missFRL5(:,winBeforeCue(end):winBeforeCue(end)+1000),2)];
    y(y<3) = NaN;
    [~,p] = ttest(y(:,1),y(:,2));
    disp(['L5 Hit vs Miss p = ' num2str(ranksum(y(:,1),y(:,2)))])
    x1 = rand(size(y,1),1)/5+1;
    x2 = rand(size(y,1),1)/5+2;
    scatter(x1,y(:,1),'b','filled')
    scatter(x2,y(:,2),'b','filled')
    xlim([.75 2.25])
    plot([x1(:)';x2(:)'], [y(:,1)';y(:,2)'], 'k-')
    title(['\color{red}L2/3 ' ' \color{blue}L5'])
    set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
    ylabel('Firing Rate (Hz)')
    
    
    
    y =  [mean(mihitFRL23(:,winBeforeCue(end):winBeforeCue(end)+1000),2),mean(mifaFRL23(:,winBeforeCue(end):winBeforeCue(end)+1000),2)];
    disp(['L23 Hit vs FA p = ' num2str(ranksum(y(:,1),y(:,2)))])
    subplot(122),hold on
    x1 = rand(size(y,1),1)/5+1;
    x2 = rand(size(y,1),1)/5+2;
    scatter(x1,y(:,1),'r','filled')
    scatter(x2,y(:,2),'r','filled')
    xlim([.75 2.25])
    plot([x1(:)';x2(:)'], [y(:,1)';y(:,2)'], 'k-')
    y =  [mean(mihitFRL5(:,winBeforeCue(end):winBeforeCue(end)+1000),2),mean(mifaFRL5(:,winBeforeCue(end):winBeforeCue(end)+1000),2)];
    y(y<3) = NaN;
    disp(['L5 Hit vs FA p = ' num2str(ranksum(y(:,1),y(:,2)))])
    x1 = rand(size(y,1),1)/5+1;
    x2 = rand(size(y,1),1)/5+2;
    scatter(x1,y(:,1),'b','filled')
    scatter(x2,y(:,2),'b','filled')
    xlim([.75 2.25])
    plot([x1(:)';x2(:)'], [y(:,1)';y(:,2)'], 'k-')
    title(['\color{red}L2/3 ' ' \color{blue}L5'])
    set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
    
    %
    figure,subplot(221),imagesc(zscore(hitFR(hitspkTuning.stimResponsive,:),[],2)),colormap(hot),caxis([0 3])
    subplot(222),imagesc(zscore(hitFR(~hitspkTuning.stimResponsive,:),[],2)),colormap(hot),caxis([0 3])
    subplot(223),imagesc(zscore(missFR(hitspkTuning.stimResponsive,:),[],2)),colormap(hot),caxis([0 3])
    subplot(224),imagesc(zscore(missFR(~hitspkTuning.stimResponsive,:),[],2)),colormap(hot),caxis([0 3])
    figure,subplot(121),
    data1 = sum(hitspkTuning.stimResponsive)/length(hitspkTuning.stimResponsive);
    data2 = 1-data1;
    pie([data1,data2])
    title('Hit')
    data1 = sum(missspkTuning.stimResponsive)/length(missspkTuning.stimResponsive);
    data2 = 1-data1;
    subplot(122),pie([data1,data2])
    title('Miss')
    legend('Responsive','Unresponsive')
    figure,
    data1 = sum(hitspkTuning.stimResponsive==1 & missspkTuning.stimResponsive==1)/length(hitspkTuning.stimResponsive(hitspkTuning.stimResponsive==1));
    data2 = 1-data1;
    pie([data1,data2])
    legend('Maintained Response','Became Unresponsive')
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

function output = spkResponse(dataTemp,winBeforeCue,winAfterCue,winRT,winRW)
for n = 1:size(dataTemp,1)
    tempFR = dataTemp(n,:);
    try
        output.stimResponsive(n) = mean(tempFR(winAfterCue(1):ceil(mean(winRT))))>(1.5*mean(tempFR(winBeforeCue))); %% responsive from stim start until reaction time
        output.supressResponsive(n) = 1.5*mean(tempFR(winAfterCue(1):ceil(mean(winRT))))<(mean(tempFR(winBeforeCue))); %% responsive from stim start until reaction time
        output.motorResponsive(n) = mean(tempFR(winAfterCue(1):ceil(nanmean(winRW))))>(1.5*mean(tempFR(ceil(nanmean(winRW)):winAfterCue(end)))); % Time from RT to reward
        output.rewardResponsive(n) = 1.5*mean(tempFR(winAfterCue(1):ceil(nanmean(winRW))))<(mean(tempFR(ceil(nanmean(winRW)):winAfterCue(end))));
    catch
        continue
    end
end
end
