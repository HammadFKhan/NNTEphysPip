temperature = data.analogChannels(1,:);
temperature = (temperature-1.25)/0.005;
IntanBehaviour.temperature = resample(temperature,parameters.Fs,data.targetedFs);
IntanBehaviour.temperature =  IntanBehaviour.temperature-IntanBehaviour.temperature(10000);
clear temperature
%%
for n = 1:IntanBehaviour.nCueHit
    hitTemp(n,1) = IntanBehaviour.temperature(IntanBehaviour.cueHitTrace(n).LFPIndex(1));
end
for n = 1:IntanBehaviour.nCueMiss
    missTemp(n,1) = IntanBehaviour.temperature(IntanBehaviour.cueMissTrace(n).LFPIndex(1));
end
for n = 1:length(IntanBehaviour.missTrace)
    FATemp(n,1) = IntanBehaviour.temperature(IntanBehaviour.missTrace(n).LFPIndex(1));
end
%%
hitRate = arrayfun(@(x) x.LFPtime(1), IntanBehaviour.cueHitTrace);
binningAmount = ceil(hitRate(end)/60);
idxh = discretize(hitRate,binningAmount);
idxc1 = zeros(1,binningAmount);
for n = 1:max(idxh)
    idxc1(n) = nanmean(hitTemp(idxh==n));
end
hitRate = histcounts(hitRate,binningAmount); % 29 should be the binning amount
%figure,plot(hitRate);
%
missRate = arrayfun(@(x) x.LFPtime(1), IntanBehaviour.cueMissTrace);
binningAmount = ceil(missRate(end)/60);
idxm = discretize(missRate,binningAmount);
idxc2 = zeros(1,binningAmount);
for n = 1:max(idxm)
    idxc2(n) = nanmean(missTemp(idxm==n));
end
missRate = histcounts(missRate,binningAmount); % 29 should be the binning amount

FARate = arrayfun(@(x) x.LFPtime(1), IntanBehaviour.missTrace);
binningAmount = ceil(FARate(end)/60);
idxf = discretize(FARate,binningAmount);
idxc3 = zeros(1,binningAmount);
for n = 1:max(idxf)
    idxc3(n) = nanmean(FATemp(idxf==n));
end
FARate = histcounts(FARate,binningAmount); % 29 should be the binning amount
idxc = nanmean([idxc1;idxc2;idxc3]);
%%
load myMap
hitRatei = interp1(1:length(hitRate),hitRate,1:0.01:length(hitRate));
idxci = interp1(1:length(idxc),idxc,1:0.01:length(idxc));
missRatei = interp1(1:length(missRate),missRate,1:0.01:length(missRate));
FARatei = interp1(1:length(FARate),FARate,1:0.01:length(FARate));
figure,
h = cline(1:length(hitRatei), smoothdata(hitRatei,'movmean',50), [], idxci);set( h, 'linestyle', '-', 'linewidth', 2  ),hold on
caxis([-20 0])
h = cline(1:length(missRatei), missRatei, [], idxci);set( h, 'linestyle', '-', 'linewidth', 2  ),hold on
caxis([-20 0])
h = cline(1:length(FARatei), FARatei, [], idxci);set( h, 'linestyle', '-', 'linewidth', 2  ),hold on
caxis([-20 0])
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
colormap(flip(myMap))
%%
figure,plot(smoothdata(hitRatei,'movmean',50),'LineWidth',2),hold on
plot(smoothdata(missRatei,'movmean',50),'LineWidth',2)
plot(smoothdata(FARatei,'movmean',50),'LineWidth',2),set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
figure,h = cline((1:length(hitRatei))*0.01, idxci, [], idxci);set( h, 'linestyle', '-', 'linewidth', 2  ),hold on
caxis([-25 0]),colormap(flip(myMap)),set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
