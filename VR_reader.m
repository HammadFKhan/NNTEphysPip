%% VR Parser
clear; clc;
[pathname,VR_data] = batchLoad();

% Distance Travelled
for i = 3:size(VR_data,2)
    maxDist(i) = max(VR_data{1,i}(:,3));
    Dist{i} = VR_data{1,i}(:,3);
    AvgVel{i} = VR_data{1,i}(:,3);
end
