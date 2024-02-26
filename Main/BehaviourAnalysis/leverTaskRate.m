function output = leverTaskRate(IntanBehaviour)
hitRate = arrayfun(@(x) x.LFPtime(1), IntanBehaviour.cueHitTrace);
binningAmount = ceil(hitRate(end)/60);
hitRate = histcounts(hitRate,binningAmount); % 29 should be the binning amount
%figure,plot(hitRate);
%
missRate = arrayfun(@(x) x.LFPtime(1), IntanBehaviour.cueMissTrace);
binningAmount = ceil(missRate(end)/60);
missRate = histcounts(missRate,binningAmount); % 29 should be the binning amount

FARate = arrayfun(@(x) x.LFPtime(1), IntanBehaviour.missTrace);
binningAmount = ceil(FARate(end)/60);
FARate = histcounts(FARate,binningAmount); % 29 should be the binning amount
%%
hitRatei = interp1(1:length(hitRate),hitRate,1:0.01:length(hitRate));
missRatei = interp1(1:length(missRate),missRate,1:0.01:length(missRate));
FARatei = interp1(1:length(FARate),FARate,1:0.01:length(FARate));
figure,plot(smoothdata(hitRatei,'movmean',50),'LineWidth',2),hold on
plot(smoothdata(missRatei,'movmean',50),'LineWidth',2)
plot(smoothdata(FARatei,'movmean',50),'LineWidth',2),set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
legend('hitRate','missRate','FARate')
%% OUTPUT structs
output.hitRate = hitRate;
output.missRate = missRate;
output.FARate = FARate;
