function Spikes = singleUnitAnalysis(fpath,VR_data)
path = [fpath,'/preAutoMerge'];
%% Read in kilosort data for matlab analysis
SpikeClusters = readNPY(fullfile(path, 'spike_clusters.npy'));
SpikeSamples = readNPY(fullfile(path, 'spike_times.npy'));
%% Analysis
Spikes.SpikeClusters = SpikeClusters;
Spikes.SpikeSamples = SpikeSamples;
Spikes = clusterSort(Spikes);
% sizePlot = ceil(sqrt(size(Spikes.Clusters,2)));
count = 1;
for i = 1:size(Spikes.Clusters,2)
    if isempty(Spikes.Clusters(i).cluster)
    else
        for ii = 1:size(Spikes.Clusters(i).cluster,2)
            x = Spikes.Clusters(i).cluster;
            y = ones(1,length(Spikes.Clusters(i).cluster));
        end
    end
end
% ISI
Spikes = ISI(Spikes,0.03,8192); %Spikes, Interval, Fs
Spikes = rateMap(Spikes,VR_data); %Trial number

% Clustered Projection
% Clusterless Projection
%% Plot Rate map
spikeRate = Spikes.VR(1).spikeRate;
norm = (spikeRate-min(spikeRate,[],1))./(max(spikeRate,[],1)-min(spikeRate,[],1));
for trial = 1:length(Spikes.VR)
    figure('name',['Spike Map Trial ' num2str(trial)]),...
    subplot(2,1,1),spikeImage = spike_map(norm',(1:Spikes.VR(trial).binWinTime*Spikes.VR(trial).time(end)));
    subplot(2,1,2),bar(Spikes.VR(trial).binWinTime*Spikes.VR.Velocity(:,1),Spikes.VR.Velocity(:,2));
    ylabel('Velocity cm/s')
    yline(mean(Spikes.VR.Velocity(:,2)),'r--'); box off
    ylim([-0.2 6])
end
%% Velocity Triggered Analysis
Velocity = Spikes.VR.Velocity(:,2);
velocityTrig = 1; %Triggered at 1 cm
loc = find(abs(Velocity)>=velocityTrig);

% Make time window around trigger

spikeRateTrig = Spikes.VR(1).spikeRate(loc,:); 
baselineSpike = mean(Spikes.VR(1).spikeRate,1);
figure
for i = 1:length(loc)
    if (loc(i)<=2 || loc(i)>=size(spikeRate,1)-2)
        spikeRatewin{i} = [];
    else
        spikeRatewin{i} = Spikes.VR(1).spikeRate(loc(i)-2:loc(i)+2,:)-baselineSpike;
    end
end
spikeRateall = (horzcat(spikeRatewin{:}));
plot(mean(spikeRateall,2));
normTrig = (spikeRateTrig-min(spikeRateTrig,[],1))./(max(spikeRateTrig,[],1)-min(spikeRateTrig,[],1));

% Plot
figure('name',['Spike Map Trial Velocity Triggered']),...
    subplot(2,1,1),spikeImage = spike_map(normTrig',(1:size(loc,1))); caxis([0 1]);
subplot(2,1,2),bar(1:size(loc,1),Spikes.VR.Velocity(loc,2));
ylabel('Velocity cm/s')
yline(velocityTrig,'r--'); box off
ylim([-0.2 6])
