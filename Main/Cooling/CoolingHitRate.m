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
hitRatei = interp1(1:length(hitRate),hitRate,1:0.01:length(hitRate));
idxci = interp1(1:length(idxc),idxc,1:0.01:length(idxc));
missRatei = interp1(1:length(missRate),missRate,1:0.01:length(missRate));
FARatei = interp1(1:length(FARate),FARate,1:0.01:length(FARate));
figure,
h = cline(1:length(hitRatei), smoothdata(hitRatei,'movmean',50), [], idxci);set( h, 'linestyle', '-', 'linewidth', 2  ),hold on
caxis([-25 0])
h = cline(1:length(missRatei), missRatei, [], idxci);set( h, 'linestyle', '-', 'linewidth', 2  ),hold on
caxis([-25 0])
h = cline(1:length(FARatei), FARatei, [], idxci);set( h, 'linestyle', '-', 'linewidth', 2  ),hold on
caxis([-25 0])
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
colormap(flip(myMap))
%%
figure,plot(smoothdata(hitRatei,'movmean',50),'LineWidth',2),hold on
plot(smoothdata(missRatei,'movmean',50),'LineWidth',2)
plot(smoothdata(FARatei,'movmean',50),'LineWidth',2),set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
figure,h = cline((1:length(hitRatei))*0.01, idxci, [], idxci);set( h, 'linestyle', '-', 'linewidth', 2  ),hold on
caxis([-25 0]),colormap(flip(myMap)),set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
