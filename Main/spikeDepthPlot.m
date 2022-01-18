function Spikes = spikeDepthPlot(Spikes,templateDepths)
temp = templateDepths(1:length(Spikes.Clusters))*25;
% find max cluster size
[depth,idx] = sort(temp);

Velocity = Spikes.VR.Velocity(:,2);
velocityTrig = 1; %Triggered at 1 cm
loc = find(abs(Velocity)>=velocityTrig);
sortedSpikeRate = Spikes.VR.spikeRate(:,idx)';

figure
minmax = (sortedSpikeRate-min(sortedSpikeRate,[],1))./(max(sortedSpikeRate,[],1)-min(sortedSpikeRate,[],1));
imagesc(1:100,depth,smoothdata(minmax,'gaussian',10)),colormap(jet),colorbar

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
ylim([0 60])
% Save output
Spikes.Depth.depth = depth;
Spikes.Depth.clusterID = idx;
Spikes.Depth.sortedSpikeRate = sortedSpikeRate;
Spikes.Depth.spikeRaster = spikeRaster;
end
