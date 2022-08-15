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
findL23 = find(depth<=300);
findL4 = find(depth>300 & depth<=400);
findL5 = find(depth>400);
L23SR = ssortedSpikeRate(findL23,:);
L4SR = ssortedSpikeRate(findL4,:);
L5SR = ssortedSpikeRate(findL5,:);

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
figure,plot(spikeRaster(:,1000:4000),'.','color','k'); 
% Save output
Spikes.Depth.depth = depth;
Spikes.Depth.clusterID = idx;
Spikes.Depth.sortedSpikeRate = sortedSpikeRate;
Spikes.Depth.spikeRaster = spikeRaster;
Spikes.spikeRate.L23SR = L23SR;
Spikes.spikeRate.L4SR = L4SR;
Spikes.spikeRate.L5SR = L5SR;
Spikes.spikeRate.L5ABSR = [L4SR;L5SR];
end
