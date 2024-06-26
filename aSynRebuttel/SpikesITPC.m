%% Loop through files
addpath(genpath('Main'));
pathname = 'Y:\Hammad\Ephys\PFFProject\RebuttelITPC'
%pathname = 'Y:\Hammad\Ephys\PFFProject\SpikeRateSTRSayan\Mon\'
filesnames = dir(fullfile(pathname,'*.mat'));
for n = 1:length(filesnames)
    load(fullfile(filesnames(n).folder,filesnames(n).name))
    if exist('itpcstats','var')
        Mice(n).ITPC = itpcstats;
    else
        Mice(n).ITPC = [];
    end
    if exist('Spikes','var')
        Mice(n).Spikes = Spikes;
    else
        Mice(n).Spikes = [];
    end
end
%%
L23Spikes = {};
L5Spikes = {};
for n = 1:length(Mice)
    L23Spikes{n} = mean(Mice(n).Spikes.spikeRate.L23SR,2);
    L5Spikes{n} = mean(Mice(n).Spikes.spikeRate.L5ABSR,2);
end
L23SpkTot = vertcat(L23Spikes{:});
L5SpkTot = vertcat(L5Spikes{:});
%%
temp = zeros(max([length(L23SpkTot) length(L5SpkTot)]),2);
temp(1:length(L23SpkTot),1) = L23SpkTot;
temp(1:length(L5SpkTot),2) = L5SpkTot;
figure(1)
clf
customBoxplot(temp(1:2:end,:))
box off,set(gca,'TickDir','out','fontsize',16),ylim([0 40])
%%
%load('\\nas01.itap.purdue.edu\puhome\desktop\ExcelData\SpikeRateCTXAnova.mat')
temp = spikeDat(:,1:6);
temp = temp(:);
temp(isnan(temp)) = [];
group = "";
mouse = [];
layer = [];
mouseNum = 5;
for n = 1:6
    group = [group;repmat(num2str(n-1),length(spikeDat(~isnan(spikeDat(:,n)),1)),1)];
    mLen = floor(length(spikeDat(~isnan(spikeDat(:,n)),1))/mouseNum);
    for nMouse = 1:mouseNum-1
        mouse = [mouse;repmat(num2str(nMouse),mLen,1)];
        %         layer = [layer;repmat('23',floor(mLen/2),1)];
        %         layer = [layer;repmat('l5',mLen-floor(mLen/2),1)];
    end
    blah = length(spikeDat(~isnan(spikeDat(:,n)),1))-(mouseNum-1)*mLen;
    mouse = [mouse;repmat(num2str(nMouse+1),blah,1)];
    %     layer = [layer;repmat('23',floor(blah/2),1)];
    %     layer = [layer;repmat('l5',blah-floor(blah/2),1)];
end

group(1,:) = [];
[p,t,stats] = anovan(temp,{group mouse},'model','interaction','varnames',{'group','mouse'});
[results,~,~,gnames] = multcompare(stats,"Dimension",[1 1]);
%% stat by mouse
for n = 1:6
    for nMouse = 1:mouseNum
        mouseTot(nMouse,n) = mean(temp(strcmp(group,num2str(n-1))& ismember(mouse,num2str(nMouse))));
    end
end
[p,t,stats] = anova2(mouseTot2);
%% ITPC
for n = 1:length(Mice)
    betaDepth{n} = Mice(n).ITPC.betaDepth;
end
figure(2),clf
%dat = vertcat(betaDepth{:});
customBoxplot(betaITPCtot)
box off,set(gca,'TickDir','out','fontsize',16),ylim([0 0.3])
%%
load('Y:\Hammad\Ephys\PFFProject\Stats\SpikeITPCRebuttel')
% itpcByMouse = cellfun(@mean,betaDepth,'UniformOutput',false);
% itpcByMouse = vertcat(itpcByMouse{:});
anova1(itpcByMouse(:,[1 3]))% Layer 5
anova1(itpcByMouse(:,[2 4]))% Layer 2/3
%% ITPC CTX and STR Fig
load('Y:\Hammad\Ephys\PFFProject\Stats\SpikeITPCCTXAnova')
dat = betaITPCCTXL23new;
mouseBin = 1:24:size(dat,1);
mouse = [];
for n = 1:length(mouseBin)
    if n == length(mouseBin)
        mouse(n,:) = mean(dat(mouseBin(n):end,:));
    else
        mouse(n,:) = mean(dat(mouseBin(n):mouseBin(n+1),:));
    end
end
%
%dat = betaITPCCTXL5new;
%dat = [dat;dat+0.05.*rand(24,6);dat+0.05.*rand(24,6),dat+0.05.*rand(24,6)];
[p,t,stats] = anova2(dat(1:72,:),24);
figure
multcompare(stats,'Estimate','column');
%%
[p,t,stats] = anova2(mouse,1);
figure
multcompare(stats,'Estimate','column');
%% W2W12
SpikesW2 = {};
SpikesW12 = {};
L5Spikes = {};
for n = 1:4
    SpikesW2{n} = [mean(Mice(n).Spikes.spikeRate.L23SR,2);mean(Mice(n).Spikes.spikeRate.L5ABSR,2)];
end
for n = 4:8
    SpikesW12{n} = [mean(Mice(n).Spikes.spikeRate.L23SR,2);mean(Mice(n).Spikes.spikeRate.L5ABSR,2)];
end
SpikesW2 = vertcat(SpikesW2{:});
SpikesW12 = vertcat(SpikesW12{:});
