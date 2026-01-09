%% Generate figure related to pooled eOPN inactivation experiments in primary motor cortex and primary motor thalamus 
% Combining eOPN data together
M1eOPN = struct();
ThalamuseOPN = struct();

files = dir(fullfile('D:\eOPNData\M1Inactivation\','*.mat'));
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    [M1eOPN(fileNum).neuralDynamics,waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);
    M1eOPN(fileNum).IntanBehaviour = IntanBehaviour;
    M1eOPN(fileNum).Spikes = Spikes;
    M1eOPN(fileNum).filename = files(fileNum).name;
end

files = dir(fullfile('D:\eOPNData\ThalamusInactivation\','*.mat'));
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    [ThalamuseOPN(fileNum).neuralDynamics,waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);
    ThalamuseOPN(fileNum).IntanBehaviour = IntanBehaviour;
    ThalamuseOPN(fileNum).filename = files(fileNum).name;
    ThalamuseOPN(fileNum).Spikes = Spikes;
end