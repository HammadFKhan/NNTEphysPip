function Spikes = clusterSort(Spikes)

for i = 1:max(Spikes.SpikeClusters)+1 %Add indices to compensate for 0 indexing
    spikeIdx = find(Spikes.SpikeClusters==i-1); %Subtract index to include 0
    Spikes.Clusters(i).clusterID = i-1; %Save cluster ID to verify match
    Spikes.Clusters(i).cluster = Spikes.SpikeSamples(spikeIdx);
end
% Deletes clusters that are counted as zero
% count = 1;
% for ii = 1:length(Spikes.Clusters)
%     if isempty(Spikes.Clusters(count).cluster)
%         Spikes.Clusters(count) = [];
%     else
%         count = count+1;
%     end
% end
