function Spikes = Spike_Detector_Single(DeltaFoverF,stddev,spikemin)


Spikes = zeros(size(DeltaFoverF,1),size(DeltaFoverF,2),size(DeltaFoverF,3));

dev = stddev*mean(std(DeltaFoverF(:,:,:),0,2));
avghx = mean(mean(DeltaFoverF(:,:,:)));
for j = 1:size(DeltaFoverF,1)
    for k = 1:size(DeltaFoverF,2)
        Spikes(j,k) = 0;
        if DeltaFoverF(j,k) >= dev + avghx && DeltaFoverF(j,k) >= spikemin
            Spikes(j,k) = 1;
        else 
        end
    end
end
