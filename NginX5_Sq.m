% clear; clc; 
% close all;
addpath(genpath('Main'));
% addpath(genpath('chronux'));
% addpath(genpath('Kilosort'));
addpath(genpath('npy-matlab'));
addpath(genpath('spikes-master'));
% IntanConcatenate legacy version
intandsFlag = 1; %make LFPs
activeElectrodes = 1:64;
chanMapFile = 'UCLA_chanmap_fixed.mat'; %UCLA Sharp
%chanMapFile = 'UCLA_chanmap_64F2.mat';
ds_filename = intanPreprocessing2(chanMapFile,intandsFlag,activeElectrodes); %IntanDs flag  %% double check file type
%% Combine intan data if needed
fpath = kilosortbinCombine();
%% Run Kilosort3 
% load only neccessary variables from memory mapped file
data = matfile(ds_filename);
fpath = data.fpath;
% Kilosort264FTestcode
Kilosort264SharpTestcode
savepath = fullfile(fpath,['loadme','.mat']);
save(savepath,'ds_filename');
clearvars -except ds_filename
%% New load me
[fname,fpath] = uigetfile();
savepath = fullfile(fpath,['loadme','.mat']);
ds_filename = fullfile(fpath,fname);
save(savepath,'ds_filename');
%% Parameters for behaviour
data = matfile(ds_filename); % ds_filename comes from loadme.mat
% check if data directory matches where the file originated; if not we note
% the new directory path

parameters.experiment = 'self'; % self - internally generated, cue - cue initiated
parameters.opto = 0; % 1 - opto ON , 0 - opto OFF
parameters.cool = 0; % No Cool 
parameters.windowBeforePull = 3; % in seconds
parameters.windowAfterPull = 2; % in seconds
parameters.windowBeforeCue = 1.5; % in seconds
parameters.windowAfterCue = 1.5; % in seconds
parameters.windowBeforeMI = 1.5; % in seconds 
parameters.windowAfterMI = 3.5; % in seconds 
parameters.delay = 0.5; %reward delay
parameters.Fs = 1000; % Eventual downsampled data
parameters.ts = 1/parameters.Fs;
parameters.IntanFs = data.targetedFs;
parameters.rows = 64;
parameters.cols = 1;

[Behaviour] = readLeverSq(parameters,data.amplifierTime);
[IntanBehaviour] = readLeverIntanSq(parameters,data.amplifierTime,data.analogChannels(1,:),data.digitalChannels,Behaviour,0);

% Calculate ITI time for trials and reward/no reward sequence
temp1 = arrayfun(@(x) x.LFPtime(1), IntanBehaviour.cueHitTrace);
temp1 = vertcat(temp1,ones(1,IntanBehaviour.nCueHit)); %  write 1 for reward given
temp2 = arrayfun(@(x) x.LFPtime(1), IntanBehaviour.cueMissTrace);
temp2 = vertcat(temp2,zeros(1,IntanBehaviour.nCueMiss)); %  write 0 for no reward given
temp = [temp1,temp2];
[~,idx] = sort(temp(1,:)); %sort by occurance
IntanBehaviour.ITI = temp(:,idx);
IntanBehaviour.parameters = parameters;
%% Plot behaviour
figure
for i=1:length(IntanBehaviour.hitTrace)
    plot(0:5000,smoothdata(IntanBehaviour.hitTrace(i).trace),'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
    try
        hitTrace(i,:) = smoothdata(IntanBehaviour.hitTrace(i).rawtrace);
    catch
        continue
    end
end
for n = 1:length(IntanBehaviour.hitTrace)
    IntanBehaviour.AvgHitTrace(n,:) = IntanBehaviour.hitTrace(n).trace;
end
IntanBehaviour.AvgHitTrace = mean(IntanBehaviour.hitTrace,1);
figure,
for n = 1:length(IntanBehaviour.missTrace)
    plot(0:5000,smoothdata(IntanBehaviour.missTrace(n).trace),'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
    missTrace(n,:) = IntanBehaviour.missTrace(n).trace;
end
IntanBehaviour.AvgMissTrace = mean(missTrace,1);
IntanBehaviour.AvgHitTrace = mean(IntanBehaviour.AvgHitTrace,1);

%% LFP probe setup for 64F and analysis
% Since there are two probes we want to seperate everything into linear
% maps for CSD and depthwise LFP analysis and then we do filtering
data = matfile(ds_filename);
%load UCLA_chanMap_64F2
load UCLA_chanmap_fixed.mat
if ~exist('lfp','var'),lfp = data.amplifierData;end
%TODO check if the field orientation during insertion is reversed (ie. probe 1 is lateral to probe 2)
probe1 = lfp(s.sorted_probe_wiring(:,5)==1,:);
probe2 = lfp(s.sorted_probe_wiring(:,5)==2,:);
chanProbe1 = s.sorted_probe_wiring(s.sorted_probe_wiring(:,5)==1,:); %needed for linear channel mapping later
chanProbe2 = s.sorted_probe_wiring(s.sorted_probe_wiring(:,5)==2,:); 
clear lfp
% LFP filter
set(0,'DefaultFigureWindowStyle','normal')
LFP.probe1= fastpreprocess_filtering(probe1,data.targetedFs);
% LFP.probe1 = bestLFP(LFP.probe1);
% LFP.probe1 = bandFilter(LFP.probe1,'depth'); % Extract LFPs based on 'depth' or 'single'
if ~isempty(probe2)
    LFP.probe2 = fastpreprocess_filtering(probe2,data.targetedFs);
end
%% Spikes analysis
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
%%% Calculate trial PSTH for lever
Spikes = leverPSTH(Spikes,IntanBehaviour);
%%% save spike output data to load into gui
savepath = fullfile(path,['spks4sorting','.mat']);
path = [fpath,'/kilosort3/' mergename];
save(savepath,'Spikes','-v7.3')
%ManualSpikeCurateGUI
%%% Basic spike analysis
% z-score spike rates
if exist('parameters','var')
    IntanBehaviour.parameters = parameters;
end
if exist('goodSpkComponents','var')
    Spikes.goodSpkComponents = unique(goodSpkComponents);
else 
    Spikes.goodSpkComponents = 1:length(Spikes.Clusters);
end
Spikes = rejectSpikes(Spikes,0.25,0.25,IntanBehaviour.parameters); % Reject spikes here for further analysis
[Spikes] = sortSpkLever(Spikes,IntanBehaviour);
[fpath,name,exts] = fileparts(ds_filename);
sessionName = [fpath,'/','Spikes.mat'];
% save(sessionName,"IntanBehaviour","fpath","parameters","-v7.3");
save(sessionName,"Spikes","IntanBehaviour","fpath","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
disp('Saved!')
