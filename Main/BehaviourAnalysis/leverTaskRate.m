function output = leverTaskRate(IntanBehaviour,behaviourOnlyFlag)
%% Calculate time-dependant behaviour performance
% behaviourOnlyFlag allows for quick analysis using only the .csv files
% instead of referencing to experiments. This is fine given that we are
% only looking at minute-minute transience
if behaviourOnlyFlag
    hitRate =  arrayfun(@(x) x.t0, IntanBehaviour.cueHitTrace)';
else
    hitRate = arrayfun(@(x) x.LFPtime(1), IntanBehaviour.cueHitTrace);
end
binningAmount = 0:60:ceil(IntanBehaviour.time(end)+60);
hitRate = histcounts(hitRate,binningAmount); % 29 should be the binning amount
%figure,plot(hitRate);
%
try
    if behaviourOnlyFlag
        missRate = arrayfun(@(x) x.t0, IntanBehaviour.cueMissTrace)';
    else
        missRate = arrayfun(@(x) x.LFPtime(1), IntanBehaviour.cueMissTrace);
    end
    
    %binningAmount = ceil(missRate(end)/60);
    missRate = histcounts(missRate,binningAmount); % 29 should be the binning amount
catch
    missRate = ones(1,30);
end

try
    
    if behaviourOnlyFlag
        FARate = arrayfun(@(x) x.t0, IntanBehaviour.missTrace)';
    else
        FARate = arrayfun(@(x) x.LFPtime(1), IntanBehaviour.missTrace);
    end
    %binningAmount = ceil(FARate(end)/60);
    FARate = histcounts(FARate,binningAmount); % 29 should be the binning amount
catch
    FARate = ones(1,30);
end
% hitRatei = interp1(1:length(hitRate),hitRate,1:0.01:length(hitRate));
% missRatei = interp1(1:length(missRate),missRate,1:0.01:length(missRate));
% FARatei = interp1(1:length(FARate),FARate,1:0.01:length(FARate));
% figure,plot(smoothdata(hitRatei,'movmean',50),'LineWidth',2),hold on
% plot(smoothdata(missRatei,'movmean',50),'LineWidth',2)
% plot(smoothdata(FARatei,'movmean',50),'LineWidth',2),set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
% legend('hitRate','missRate','FARate')
%% OUTPUT structs
output.hitRate = hitRate;
output.missRate = missRate;
output.FARate = FARate;
