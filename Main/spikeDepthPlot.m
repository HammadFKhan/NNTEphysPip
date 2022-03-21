function Spikes = spikeDepthPlot(Spikes,templateDepths)
temp = templateDepths(1:length(Spikes.Clusters))*25;
% find max cluster size
[depth,idx] = sort(temp);
% Removes nan values from depth
depth(isnan(depth)) = [];
Velocity = Spikes.VR.Velocity(:,2);
velocityTrig = 1; %Triggered at 1 cm
loc = find(abs(Velocity)>=velocityTrig);
sortedSpikeRate = Spikes.VR.spikeRate(:,idx)';

figure
minmax = (sortedSpikeRate-min(sortedSpikeRate,[],1))./(max(sortedSpikeRate,[],1)-min(sortedSpikeRate,[],1));
ssortedSpikeRate = smoothdata(sortedSpikeRate,'gaussian',10);
imagesc(1:100,depth,smoothdata(sortedSpikeRate,'gaussian',10)),colormap(jet),colorbar
% Spike rate by depth (indescriminate of running or no running)
findL23 = find(depth<250);
findL4 = find(depth>250 & depth<400);
findL5 = find(depth>400);
L23avgSR = mean(ssortedSpikeRate(findL23,:),2);
L4avgSR = mean(ssortedSpikeRate(findL4,:),2);
L5avgSR = mean(ssortedSpikeRate(findL5,:),2);

for i  = 1:length(Spikes.Clusters)
    clustSize(i) = size(Spikes.Clusters(i).spikeTime,2);
end

spikeRaster = zeros(length(clustSize),max(clustSize));
for i = 1:length(Spikes.Clusters)
    spikeCell = Spikes.Clusters(i).spikeTime;
    pad = abs(length(spikeCell)-max(clustSize));
    spikeRaster(i,:) = [spikeCell zeros(1,pad)];
end
spikeRaster(spikeRaster==0) = nan;
figure,plot(spikeRaster,'.','color','k'),view([90 -90]); 
ylim([5 10])
% Save output
Spikes.Depth.depth = depth;
Spikes.Depth.clusterID = idx;
Spikes.Depth.sortedSpikeRate = sortedSpikeRate;
Spikes.Depth.spikeRaster = spikeRaster;
Spikes.spikeRate.L23avgSR = L23avgSR;
Spikes.spikeRate.L4avgSR = L4avgSR;
Spikes.spikeRate.L5avgSR = L5avgSR;
end
