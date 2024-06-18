function [Spikes] = sortSpkLever(Spikes,Behaviour)
showplot = 1;
%%% zscore normalize
% hits
temp = smoothdata(Spikes.PSTH.hit.spkRates,2,'gaussian',40); %template plus smoothing
[hitnormSpk,hittimIdx,hitspkIdx] = spknorm(temp);
% misses
temp = smoothdata(Spikes.PSTH.miss.spkRates,2,'gaussian',40); %template plus smoothing
[missnormSpk,misstimIdx,missspkIdx] = spknorm(temp);
% MI hits
temp = smoothdata(Spikes.PSTH.MIHit.spkRates,2,'gaussian',40); %template plus smoothing
[MIhitnormSpk,MIhittimIdx,MIhitspkIdx] = spknorm(temp);
% MI false alarms
temp = smoothdata(Spikes.PSTH.MIFA.spkRates,2,'gaussian',40); %template plus smoothing
[MIFAnormSpk,MIFAtimIdx,MIFAspkIdx] = spknorm(temp);

%%% jitter response and stability
hitPath = hittimIdx(hitspkIdx);
missPath = misstimIdx(missspkIdx);
MIhitPath = MIhittimIdx(MIhitspkIdx);
MIFAPath = MIFAtimIdx(MIFAspkIdx);
hitJitter = diff(hitPath);
missJitter = diff(missPath);
%%% Spike Depth Assignment TODO: Gamma GED spike coherence
l23idx = arrayfun(@(x) x.spikeDepth<400,Spikes.Clusters,'UniformOutput',false); % L23
l5idx = arrayfun(@(x) x.spikeDepth>=400,Spikes.Clusters,'UniformOutput',false); % L5
l23idx = find(cell2mat(l23idx)==1);
l5idx = find(cell2mat(l5idx)==1);
% hitL23spkIdx = findLayer(Spikes.PSTH.hit.spkIdx,l23idx);
% missL23spkIdx = findLayer(Spikes.PSTH.miss.spkIdx,l23idx);
% hitL5spkIdx = findLayer(Spikes.PSTH.hit.spkIdx,l5idx);
% missL5spkIdx = findLayer(Spikes.PSTH.miss.spkIdx,l5idx);
% 
% MIhitL23spkIdx = findLayer(Spikes.PSTH.MIHit.spkIdx,l23idx);
% MIFAL23spkIdx = findLayer(Spikes.PSTH.MIFA.spkIdx,l23idx);
% MIhitL5spkIdx = findLayer(Spikes.PSTH.MIHit.spkIdx,l5idx);
% MIFAL5spkIdx = findLayer(Spikes.PSTH.MIFA.spkIdx,l5idx);
% MIhitL23Path = MIhittimIdx(MIhitL23spkIdx);
% MIhitL5Path = MIhittimIdx(MIhitL5spkIdx);
% MIFAL23Path = MIFAtimIdx(MIFAL23spkIdx);
% MIFAL5Path = MIFAtimIdx(MIFAL5spkIdx);

%%% plot it out
% cell map
if showplot
    %%% All neurons
    figure,subplot(121),imagesc(-Behaviour.parameters.windowBeforeCue*1000:Behaviour.parameters.windowAfterCue*1000,...
        1:size(hitnormSpk,1),hitnormSpk(hitspkIdx,:)),hold on
    colormap(flip(gray))
    colorbar
    set(gca,'fontsize',16)
    caxis([0.0 2.56])
    RT = mean(Behaviour.reactionTime)*1000;
    xline(0);
    xline(RT);
    subplot(121),plot((hitPath)-Behaviour.parameters.windowBeforeCue*1000,1:size(hitnormSpk,1),'r','LineWidth',1)
    title('Cue Hit')
    subplot(122),imagesc(-Behaviour.parameters.windowBeforeCue*1000:Behaviour.parameters.windowAfterCue*1000,...
        1:size(missnormSpk,1),missnormSpk(missspkIdx,:)),hold on
    colormap(flip(gray))
    colorbar
    set(gca,'fontsize',16)
    caxis([0.0 2.56])
    RT = mean(Behaviour.reactionTime)*1000;
    xline(0);
    xline(RT);
    subplot(122),plot((missPath)-Behaviour.parameters.windowBeforeCue*1000,1:size(missnormSpk,1),'r','LineWidth',1)
    title('Cue Miss')
    %%% Motion Initiated
    
     figure,subplot(121),imagesc(-Behaviour.parameters.windowBeforeCue*1000:Behaviour.parameters.windowAfterCue*1000,...
        1:size(MIhitnormSpk,1),MIhitnormSpk(MIhitspkIdx,:)),hold on
    colormap(flip(gray))
    colorbar
    set(gca,'fontsize',16)
    caxis([0.0 2.56])
    xline(0);
    subplot(121),plot(MIhitPath-Behaviour.parameters.windowBeforeCue*1000,1:size(MIhitnormSpk(MIhitspkIdx),1),'r','LineWidth',1)
    title('MI Hit')
    
    subplot(122),imagesc(-Behaviour.parameters.windowBeforeCue*1000:Behaviour.parameters.windowAfterCue*1000,...
        1:size(MIFAnormSpk,1),MIFAnormSpk(MIFAspkIdx,:)),hold on
    colormap(flip(gray))
    colorbar
    set(gca,'fontsize',16)
    caxis([0.0 2.56])
    
%     RT = mean(Behaviour.reactionTime)*1000;
    xline(0);
    xline(RT);
    subplot(122),plot(MIFAPath-Behaviour.parameters.windowBeforeCue*1000,1:size(MIFAnormSpk(MIFAspkIdx),1),'r','LineWidth',1)
    title('MI False Alarm')
    
    
    
    % Jitter response
    figure,bar(1,mean(hitJitter),'FaceColor',[188/255 190/255 192/255]),hold on
    errorbar(1,mean(hitJitter),std(hitJitter)/sqrt(length(hitJitter)),'k')
    scatter(1*ones(length(hitJitter(hitJitter~=0)),1),hitJitter(hitJitter~=0),'b','filled','jitter','on','jitterAmount',0.1)
    
    bar(2,mean(missJitter),'FaceColor',[188/255 190/255 192/255]),hold on
    errorbar(2,mean(missJitter),std(missJitter)/sqrt(length(missJitter)),'k')
    scatter(2*ones(length(missJitter(missJitter~=0)),1),missJitter(missJitter~=0),'b','filled','jitter','on','jitterAmount',0.1)
    box off
    set(gca,'tickdir','out')
    set(gca,'fontsize',16)
    ylabel('Assembly Jitter')
    ylim([0 100])
end
%%% output
Spikes.PSTH.hit.normSpk = hitnormSpk;
Spikes.PSTH.hit.spkIdx = hitspkIdx;
Spikes.PSTH.miss.normSpk = missnormSpk;
Spikes.PSTH.miss.spkIdx = missspkIdx;
% Spikes.PSTH.hit.l23spkidx = hitl23spkidx;
% Spikes.PSTH.miss.l23spkidx = missl23spkidx;
% Spikes.PSTH.hit.l5spkidx = hitl5spkidx;
% Spikes.PSTH.miss.l5spkidx = missl5spkidx;
Spikes.PSTH.MIHit.normSpk = MIhitnormSpk;
Spikes.PSTH.MIHit.spkIdx = MIhitspkIdx;
Spikes.PSTH.MIFA.normSpk = MIFAnormSpk;
Spikes.PSTH.MIFA.spkIdx = MIFAspkIdx;
end

function output = findLayer(spkIdx,idx)
output = ismember(spkIdx,idx);
output = spkIdx(output);
end
function [normSpikeRate,idx,idxc] = spknorm(temp)
[nanIdx,~,~] = find(~isnan(temp));
nanIdx = unique(nanIdx);
normSpikeRate = zscore(temp(nanIdx,:),0,2);
idx = zeros(size(normSpikeRate,1),1);
for n = 1:length(idx)
    [~,idx(n)] = max(normSpikeRate(n,:));
end
[~,idxc] = sort(idx);
end