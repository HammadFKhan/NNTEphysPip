%% Kilosort data preprocessing
datI = int16(amplifier_data);
kilosortOut = uigetdir();
[SUCCESS,~,~] = mkdir(kilosortOut,filename(1:end-4));
newDIR = dir(strcat(kilosortOut,'\',filename(1:end-4)));
newDIR = newDIR.folder;
fid = fopen(strcat(newDIR,'\',filename(1:end-4)),'w');
fwrite(fid,datI,'int16');
fclose(fid);
%% Kilosort Analysis

path = uigetdir();
SpikeClusters = readNPY(fullfile(path, 'spike_clusters.npy'));
SpikeSamples = readNPY(fullfile(path, 'spike_times.npy'));

Spikes.SpikeClusters = SpikeClusters;
Spikes.SpikeSamples = SpikeSamples;
Spikes = clusterSort(Spikes);
sizePlot = ceil(sqrt(size(Spikes.Clusters,2)));
for i = 1:size(Spikes.Clusters,2)
    for ii = 1:size(Spikes.Clusters(i).cluster,2)
        x = Spikes.Clusters(i).cluster;
        y = ones(1,length(Spikes.Clusters(i).cluster));
    end
    subplot(sizePlot,sizePlot,i),scatter(x,y,3), axis tight, box off;
end

%% ISI
Spikes = ISI(Spikes,0.005);
Spikes = rateMap(Spikes,VR_data);
%%
figure,spikeImage = spike_map(Spikes.VR.spikeCount',(1:67)*2);
