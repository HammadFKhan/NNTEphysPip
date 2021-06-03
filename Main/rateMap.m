function Spikes = rateMap(Spikes,VR_data)

for trial = 1:length(VR_data.Position)
    position = VR_data.Position{trial};
    edgesPos = 0:1:200; %Bin width on the track
    VRposition = discretize(position,edgesPos);
    VRtime = VR_data.Time{trial};
    for i = 1:length(Spikes.Clusters)
        for ii = 1:VRposition(end)
            binFind = find(VRposition==ii);
            map(ii,1) = min(binFind);
            map(ii,2) = max(binFind);
            
            if trial-1 == 0
                mapTime(ii,1) = VRtime(map(ii,1));
                mapTime(ii,2) = VRtime(map(ii,2));
            else
                mapTime(ii,1) = VRtime(map(ii,1)) + max(VR_data.Time{1,trial-1}); %adds the previous trial if ther was one
                mapTime(ii,2) = VRtime(map(ii,2)) + max(VR_data.Time{1,trial-1});
            end
            spikeCount(ii,1) = length(find(mapTime(ii,1)<Spikes.Clusters(i).spikeTime...
                & Spikes.Clusters(i).spikeTime<mapTime(ii,2)));
            
        end
        spikeRate(:,i) = spikeCount./(mapTime(:,2)-mapTime(:,1));
        Spikes.VR(trial).spikeCount(:,i) = spikeCount;
    end
    
    Spikes.VR(trial).time = VR_data.Time;
    Spikes.VR(trial).mapTime = mapTime;
    Spikes.VR(trial).position = VRposition;
    Spikes.VR(trial).spikeRate = spikeRate;
    clear mapTime spikeRate VRposition spikeCount
end