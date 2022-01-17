function spikeDepthPlot(Spikes,templateDepths)
temp = templateDepths(1:length(Spikes.Clusters))*25;
% find max cluster size
[idx,loc] = sort(temp);
norm = (Spikes.VR.spikeRate-min(spikeRate,[],1))./(max(Spikes.VR.spikeRate,[],1)-min(Spikes.VR.spikeRate,[],1));
imagesc(norm(:,loc)),colormap(jet)
figure,
for i = 1:5
    disp(i)
    plot(Spikes.Clusters(i).spikeTime,i,'.'),hold on;
end
