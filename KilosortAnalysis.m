%% Kilosort Analysis
addpath(genpath('main'));
path = uigetdir();
SpikeClusters = readNPY(fullfile(path, 'spike_clusters.npy'));
SpikeSamples = readNPY(fullfile(path, 'spike_times.npy'));
%% Analysis
Spikes.SpikeClusters = SpikeClusters;
Spikes.SpikeSamples = SpikeSamples;
Spikes = clusterSort(Spikes);
sizePlot = ceil(sqrt(size(Spikes.Clusters,2)));
for i = 1:size(Spikes.Clusters,2)
    for ii = 1:size(Spikes.Clusters(i).cluster,2)
        x = Spikes.Clusters(i).cluster;
        y = ones(1,length(Spikes.Clusters(i).cluster));
    end
    subplot(sizePlot,sizePlot,i),scatter(x,y,3), axis tight, box off;
end
% ISI
Spikes = ISI(Spikes,0.005);
Spikes = rateMap(Spikes,VR_data);
% Clustered Projection
% CLusterless Projection
%% Plot All
figure('name','Spike Map'),spikeImage = spike_map(Spikes.VR.spikeCount',(1:67)*2);
