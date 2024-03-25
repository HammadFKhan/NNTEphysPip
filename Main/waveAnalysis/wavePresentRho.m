function [waveDir,waveRho,rhoMax,rhoDiff] = wavePresentRho(Waves)
%% Plot wave Rho value per wave present scatter plot
wavePresent = arrayfun(@(x) x.wavePresent,Waves.wavesHit,'UniformOutput',false);
wavePresent = vertcat(wavePresent{:});
waveStart = arrayfun(@(x) x.waveStart,Waves.wavesHit,'UniformOutput',false);
waveStart = vertcat(waveStart{:});
%% Discribe linearity
waveRholin = arrayfun(@(x) x.rho, Waves.wavesHit, 'UniformOutput', false);
figure,
count = 1;
for n = 1:length(waveRholin)
    for nn = 1:length(waveRholin{n})
        rhoMax(count) = max(waveRholin{n}{nn});
        rhoDiff(count) = max(diff(waveRholin{n}{nn}));
        scatter(rhoMax(count),rhoDiff(count),'k'),hold on
        count = count+1;
    end
end
%%
for n = 1:size(wavePresent,1)
    dirTemp = zeros(1,3001);
    rhoTemp = zeros(1,3001);
    idx = find(waveStart(n,:)==1);
    idx2 = find(wavePresent(n,:)==1);
    temp = vertcat(Waves.wavesHit(n).waveDir);
    temp2 = vertcat(Waves.wavesHit(n).rho{:});
    dirTemp(idx) = temp(1:length(idx));
    rhoTemp(idx2) = temp2(1:length(idx2));
    waveDir(n,:) = dirTemp;
    waveRho(n,:) = rhoTemp;
end
waveDir(waveDir==0) = NaN;
waveRho(waveRho==0) = NaN;
load myMap
tIdx = randi(size(waveStart,1),60,1);
waveDirt = waveDir(tIdx,:);
waveStartt = waveStart(tIdx,:);
wavePresentt = wavePresent(tIdx,:);
waveRhot= waveRho(tIdx,:);
figure,
for n = 1:size(waveStartt,1)
    scatter(1:3001,n*waveStartt(n,:),10,waveDirt(n,:),'filled'),hold on,axis tight,colormap(myMap)
end
figure,
for n = 1:size(wavePresentt,1)
    scatter(1:3001,n*wavePresentt(n,:),10,waveRhot(n,:),'filled'),hold on,axis tight,colormap(myMap)
end
