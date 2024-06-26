load Y:\Austin\ShankData\240612_shank_eOPN_M1\loadme
load Y:\Austin\ShankData\240612_shank_eOPN_M1\IntanBehaviour

[fpath,name,exts] = fileparts(ds_filename);
data = matfile(ds_filename);
path = [fpath,'/kilosort3/'];
mergename = 'merged';
Kilosort3AutoMergeTester
path = [fpath,'/kilosort3/' mergename];
%%% Spike preprocessing (includes merging (optional) and channel info
%%% return)
% Read in kilosort data for matlab analysis
SpikeClusters = readNPY(fullfile(path, 'spike_clusters.npy'));
SpikeSamples = readNPY(fullfile(path, 'spike_times.npy'));
SpikeChannel = readNPY(fullfile(path,'channel_positions.npy'));
Spikes.SpikeClusters = SpikeClusters; 
Spikes.SpikeSamples = SpikeSamples;
Spikes = clusterSort(Spikes);
Spikes = ISI(Spikes,0.01,data.Fs,0); %Spikes, Interval, Fs
% Calculate Depth profile
%load chanMap64F2
load chanMap64Sharp
[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms, max_site] =...
    spikeTemplatePosition(data.fpath,ycoords,[]); % 'invert'
for i = 1:length(tempAmps)
    Spikes.Clusters(i).spikeDepth = templateDepths(i);
    Spikes.Clusters(i).channelDepth = max_site(i);
    Spikes.Clusters(i).spikeAmplitude = tempAmps(i);
    Spikes.Clusters(i).waveforms = waveforms(i,:);
    Spikes.Clusters(i).spikeDuration = templateDuration(i)/data.Fs*1000;
end
%%% delete empty spikes
temp = arrayfun(@(x) isempty(x.cluster), Spikes.Clusters);
Spikes.Clusters(temp) = []; 
%% Calculate trial PSTH for passive stimulation via pole

% TODO: Make a parameters structure to track experimental settings. Need
% this for correct analysis
parameters.windowBeforePole = 0.5;
parameters.windowAfterPole = 0.75;
parameters.experiment = 'passive'; % passive - pole, active - active touch
parameters.opto = 0; % 1 - opto ON , 0 - opto OFF
parameters.Fs = 1000; % Eventual downsampled data
parameters.ts = 1/parameters.Fs;

Spikes = polePSTH(Spikes,IntanBehaviour);
%%% save spike output data to load into gui
savepath = fullfile(path,['spks4sorting','.mat']);
save(savepath,'Spikes','-v7.3')
%% Basic spike analysis
% z-score spike rates
if exist('parameters','var')
    IntanBehaviour.parameters = parameters;
end
if exist('goodSpkComponents','var')
    Spikes.goodSpkComponents = unique(goodSpkComponents);
else 
    Spikes.goodSpkComponents = 1:length(Spikes.Clusters);
end
Spikes = rejectSpikespassive(Spikes,0.25,0.25,parameters); % Reject spikes here for further analysis
[fpath,name,exts] = fileparts(ds_filename);
sessionName = [fpath,'/','Spikes.mat'];
% save(sessionName,"IntanBehaviour","fpath","parameters","-v7.3");
save(sessionName,"Spikes","IntanBehaviour","fpath","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
disp('Saved!')

%% Plot some rep units
figure,hold on
time = -parameters.windowBeforePole*parameters.Fs:parameters.windowAfterPole*parameters.Fs;
for n = 3
subplot(2,1,[1]),Show_Spikes(Spikes.PSTH.pole.spks{n}),axis off
subplot(2,1,[2]),bar(time,smoothdata(Spikes.PSTH.pole.spkRates(n,:)),'FaceColor',[28/255 117/255 188/255],'EdgeColor','none')
axis tight, box off, set(gca,'TickDir','out')
set(gca,'fontsize',16)
xline(0)
end