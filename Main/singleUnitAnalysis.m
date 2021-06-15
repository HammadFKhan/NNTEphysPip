function Spikes = singleUnitAnalysis(fpath,VR_data)
path = [fpath,'/postAutoMerge'];
%% Read in kilosort data for matlab analysis
SpikeClusters = readNPY(fullfile(path, 'spike_clusters.npy'));
SpikeSamples = readNPY(fullfile(path, 'spike_times.npy'));

%% Analysis
Spikes.SpikeClusters = SpikeClusters;
Spikes.SpikeSamples = SpikeSamples;
Spikes = clusterSort(Spikes);
% sizePlot = ceil(sqrt(size(Spikes.Clusters,2)));
count = 1;
for i = 1:size(Spikes.Clusters,2)
    if isempty(Spikes.Clusters(i).cluster)
    else
        for ii = 1:size(Spikes.Clusters(i).cluster,2)
            x = Spikes.Clusters(i).cluster;
            y = ones(1,length(Spikes.Clusters(i).cluster));
        end
    end
end
% ISI
Spikes = ISI(Spikes,0.03);
Spikes = rateMap(Spikes,VR_data); %Trial number

% Clustered Projection
% Clusterless Projection
%% Plot Rate map
for trial = 1:length(Spikes.VR)
    figure('name',['Spike Map Trial ' num2str(trial)]),...
        spikeImage = spike_map(Spikes.VR(trial).spikeRate',(1:Spikes.VR(trial).position(end)));
%     subplot(2,1,2),plot(smoothdata(VR_data.AvgVel{1,trial},'sgolay')); axis tight, axis off;
end
