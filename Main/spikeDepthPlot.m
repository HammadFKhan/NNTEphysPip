function spikeDepthPlot(Spikes,templateDepths)
temp = templateDepths*2*25;
% find max cluster size
clusterMax = max(Spikes.Clusters.spikeTime)
figure,
for i = 1:5
    disp(i)
    plot(Spikes.Clusters(i).spikeTime,i,'.'),hold on;
end
