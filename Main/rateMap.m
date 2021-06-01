function Spikes = rateMap(Spikes,VR_Data,varargin)
if nargin<3 || strcmp(trialNumber,''), display('No Trial number was selected!') 
end

p = inputParser;
trialNumber = 3;

addParameter(p,'trialNumber',trialNumber,@isnumeric);
parse(p,Spikes,VR_Data,varargin{:});

trialNumber = p.Results.trialNumber;

position = VR_Data{1,trialNumber}(:,3);
edgesPos = 0:2:200;
VRposition = discretize(position,edgesPos);
VRtime = VR_Data{1,trialNumber}(:,1);
for i = 1:length(Spikes.Clusters)
    for ii = 1:VRposition(end)
        binFind = find(VRposition==ii);
        map(ii,1) = min(binFind);
        map(ii,2) = max(binFind);
        mapTime(ii,1) = VRtime(map(ii,1));
        mapTime(ii,2) = VRtime(map(ii,2));
        spikeCount(ii,1) = length(find(mapTime(ii,1)<Spikes.Clusters(i).spikeTime...
            & Spikes.Clusters(i).spikeTime<mapTime(ii,2)));
    end
    Spikes.VR.spikeCount(:,i) = spikeCount;
    Spikes.VR.VR_data = VR_Data;
    Spikes.VR.mapTime = mapTime;
    Spikes.VR.position = VRposition;
end