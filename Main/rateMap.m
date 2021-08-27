function Spikes = rateMap(Spikes,VR_data)

for trial = 1:length(VR_data.Position)
    position = VR_data.Position{trial};
    edgesPos = 0:1:ceil(VR_data.Position{1}(end));
    VRposition =  discretize(position,edgesPos);
    time = VR_data.Time{trial};
    edgesTime = 0:3:16600; %Bin width on the track
    VRtimebin = discretize(time,edgesTime);
    checkBin = find(diff(VRtimebin)>1);
    
    if isempty(size(checkBin,1))==0 %Checks for binning error
        for check = 1:size(checkBin,1)
            checkValue = checkBin(check)+1; %Adds 1 to the location of the skipped cell
            VRtimebin(checkValue,1) = VRtimebin(checkValue,1)-1;
        end
    end
    VRtime = VR_data.Time{trial};
    dP = diff(position);
    dP = vertcat(0,dP);
    for i = 1:length(Spikes.Clusters)
        for ii = 1:VRtimebin(end)
            binFind = find(VRtimebin==ii);
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
    % Custom moving velocity window
    for i = 1:VRtimebin(end)
        findBinT = find(VRtimebin==i);
        sumPosition = sum(dP(findBinT,1)); %#ok<FNDSB>
        Vel(i,1) = i;
        Vel(i,2) = sumPosition;
    end

    spikeRate(spikeRate>=200)=0; % Negates outliers
    spikeRate(isnan(spikeRate))=0; % If spikeRate contains NaN
    Spikes.VR(trial).time = VRtimebin;
    Spikes.VR(trial).mapTime = mapTime;
    Spikes.VR(trial).position = VRposition;
    Spikes.VR(trial).spikeRate = spikeRate;
    Spikes.VR(trial).Velocity = Vel;
    clear mapTime spikeRate VRposition spikeCount
end