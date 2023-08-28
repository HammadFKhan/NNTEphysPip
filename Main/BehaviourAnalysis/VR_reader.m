%% VR Parser
addpath(genpath('main'));
[pathname,VRdat] = batchLoad();


% Distance Travelled
count = 1;
for i = 3:size(VRdat,2)
    Dist{count} = max(VRdat{1,i}(:,3));
    Position{count} = VRdat{1,i}(:,3);
    AvgVel{count} = VRdat{1,i}(:,5);
    Time{count} = VRdat{1,i}(:,1);
    count = count+1;
end
VR_data.AvgVel = AvgVel;
VR_data.Position = Position;
VR_data.Dist = vertcat(Dist{:});
VR_data.Time = Time;

clear Time AvgVel Position Dist count