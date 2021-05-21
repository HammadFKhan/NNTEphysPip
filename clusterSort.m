function Spikes = clusterSort(Spikes)
% Parse
SpikeClusters = Spikes.SpikeClusters;
SpikeSamples = Spikes.SpikeSamples;
for i = 1:max(SpikeClusters)
    count = 1;
    for ii = 1:size(SpikeSamples,1)
        if SpikeClusters(ii,:) == i
            Spikes.Clusters(i).cluster(count) = SpikeSamples(ii,:);
            count = count+1;
        end
    end
end
