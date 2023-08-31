% clear; clc; 
% close all;
addpath(genpath('Main'));
% addpath(genpath('chronux'));
% addpath(genpath('Kilosort'));
addpath(genpath('npy-matlab'));
addpath(genpath('spikes-master'));
% IntanConcatenate legacy version
ds_filename = intanPreprocessing2;
% Run Kilosort3 
% load only neccessary variables from memory mapped file
data = matfile(ds_filename);
fpath = data.fpath;
Kilosort264FTestcode
savepath = fullfile(fpath,['loadme','.mat']);
save(savepath,'ds_filename');
clearvars -except ds_filename
%% Parameters for behaviour
parameters.experiment = 'cue'; % self - internally generated, cue - cue initiated
parameters.opto = 0; % 1 - opto ON , 0 - opto OFF
parameters.windowBeforePull = 1; % in seconds
parameters.windowAfterPull = 1; % in seconds
parameters.windowBeforeCue = 0.5; % in seconds
parameters.windowAfterCue = 1.5; % in seconds
parameters.Fs = 1000;
parameters.ts = 1/parameters.Fs;
[Behaviour] = readLever(parameters,data.amplifierTime);
[IntanBehaviour] = readLeverIntan(parameters,data.amplifierTime,data.analogChannels,data.digitalChannels,Behaviour);
%% Plot behaviour
figure
for i=1:IntanBehaviour.nCueHit
    plot(0:2000,smoothdata(IntanBehaviour.cueHitTrace(i).trace),'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
    try
        hitTrace(i,:) = smoothdata(IntanBehaviour.cueHitTrace(i).rawtrace);
    catch
        continue
    end
end
for n = 1:IntanBehaviour.nCueHit
    IntanBehaviour.AvgHitTrace(n,:) = IntanBehaviour.cueHitTrace(n).trace;
end
IntanBehaviour.AvgHitTrace = mean(IntanBehaviour.AvgHitTrace,1);
for n = 1:IntanBehaviour.nCueMiss
    IntanBehaviour.AvgMissTrace(n,:) = IntanBehaviour.cueMissTrace(n).trace;
end
IntanBehaviour.AvgMissTrace = mean(IntanBehaviour.AvgMissTrace,1);
%% LFP probe setup for 64F and analysis
% Since there are two probes we want to seperate everything into linear
% maps for CSD and depthwise LFP analysis and then we do filtering
data = matfile(ds_filename);
load UCLA_chanMap_64F
if ~exist('lfp','var'),lfp = data.amplifierData;end
%TODO check if the field orientation during insertion is reversed (ie. probe 1 is lateral to probe 2)
probe1 = lfp(s.sorted_probe_wiring(:,5)==1,:);
probe2 = lfp(s.sorted_probe_wiring(:,5)==2,:);
chanProbe1 = s.sorted_probe_wiring(s.sorted_probe_wiring(:,5)==1,:); %needed for linear channel mapping later
chanProbe2 = s.sorted_probe_wiring(s.sorted_probe_wiring(:,5)==2,:); 

%% LFP filter
set(0,'DefaultFigureWindowStyle','normal')
LFP.probe1 = fastpreprocess_filtering(probe1,data.targetedFs);
LFP.probe1 = bestLFP(LFP.probe1);
% LFP.probe1 = bandFilter(LFP.probe1,'depth'); % Extract LFPs based on 'depth' or 'single'
LFP.probe2 = fastpreprocess_filtering(probe2,data.targetedFs);
LFP.probe2 = bestLFP(LFP.probe2);
% LFP.probe2 = bandFilter(LFP.probe2,'depth'); % Extract LFPs based on 'depth' or 'single'
% Build linear channels (should be four from the 64F)
LFP.probe1.chan1 = find(chanProbe1(:,2)==300);
LFP.probe1.chan2 = find(chanProbe1(:,2)==320.1);
LFP.probe2.chan1 = find(chanProbe2(:,2)==0);
LFP.probe2.chan2 = find(chanProbe2(:,2)==20);
%% CSD and spectrogram
LFP.probe1.power = leverLFPAnalysis(LFP.probe1.LFP(LFP.probe1.chan2,:),IntanBehaviour);
%% Spikes analysis
path = [data.fpath,'/kilosort3'];
% Read in kilosort data for matlab analysis
SpikeClusters = readNPY(fullfile(path, 'spike_clusters.npy'));
SpikeSamples = readNPY(fullfile(path, 'spike_times.npy'));
Spikes.SpikeClusters = SpikeClusters; %Add one because of 0 index from python
Spikes.SpikeSamples = SpikeSamples;
Spikes = clusterSort(Spikes);
Spikes = ISI(Spikes,0.01,data.Fs,0); %Spikes, Interval, Fs
% Calculate Depth profile
load chanMap64F
[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] =...
    spikeTemplatePosition(data.fpath,ycoords);
for i = 1:length(tempAmps)
    Spikes.Clusters(i).spikeDepth = templateDepths(i);
    Spikes.Clusters(i).spikeAmplitude = tempAmps(i);
    Spikes.Clusters(i).waveforms = waveforms(i,:);
    Spikes.Clusters(i).spikeDuration = templateDuration(i)/data.Fs*1000;
end
%% Calculate trial PSTH for lever
Spikes = leverPSTH(Spikes,Behaviour);
%% Neural Trajectory Analysis using GPFA
Spikes = makeSpikeGPFA(Spikes);
[result,seqTrain] = gpfaAnalysis(Spikes.GPFA.hit.dat);
%% Some local function to make things easier 
% Hit trials
function output = leverLFPAnalysis(linearProbe,Behaviour) % LFP of linear channel and Behaviour struct
CSDoutputhit = [];waveletHit = [];waveletMiss = [];powerCWThit = [];CSDoutputmiss = []; hitLFP = [];missLFP = [];
params.tapers = [5 9];
params.Fs = 1000;
params.fpass = [4 80];
params.err = [2 0.05];
for i = 1:Behaviour.nCueHit
    timestamphit(i,:) = [Behaviour.cueHitTrace(i).LFPtime(1),Behaviour.cueHitTrace(i).LFPtime(end)]; %taken in seconds
    hitWin = floor(Behaviour.cueHitTrace(i).LFPtime*1000); %multiply by Fs;
    hitLFP(:,:,i) = linearProbe(:,hitWin);
    [powerCWThit(:,:,i), fwt] = calCWTSpectogram(mean(hitLFP(1:2,:,i),1),0:2000,1000,10,[10 40],0,0);
%     [CSDoutputhit(:,:,i)]  = CSD(hitLFP(:,:,i)/1E6,1000,50E-6);
end
for i = 1:Behaviour.nCueMiss
    timestampmiss(i,:) = [Behaviour.cueMissTrace(i).LFPtime(1),Behaviour.cueMissTrace(i).LFPtime(end)]; %taken in seconds
    missWin = floor(Behaviour.cueMissTrace(i).LFPtime*1000);
    missLFP(:,:,i) = linearProbe(:,missWin);
    [powerCWTmiss(:,:,i), fwt] = calCWTSpectogram(mean(missLFP(1:2,:,i),1),0:2000,1000,10,[10 40],0);
    %     [CSDoutputmiss(:,:,i)]  = CSD(missLFP(:,:,i)'/1E6,1000,20E-6);
end
for i = 1:size(linearProbe,1)
    tfmiss(:,:,i) = itpc(linearProbe(i,:),timestampmiss,1000,0);
    [tfhit(:,:,i),frex,pnts] = itpc(linearProbe(i,:),timestamphit,1000,0);
end
% figure(), clf
% contourf(1:pnts,frex,mean(tfhit,3),10,'linecolor','none')
% set(gca,'clim',[0 .2],'ydir','normal')
% title('ITPC')
% colormap(jet),colorbar
% figure,plot(smoothdata(squeeze(mean(tfhit,1)),'gaussian',50),'k')
% hold on,plot(smoothdata(squeeze(mean(tfmiss,1)),'gaussian',50),'r')
[Shit,fhit,Serrhit]=mtspectrumc(squeeze(mean(hitLFP,1)),params);
[Smiss,~,Serrmiss] = mtspectrumc(squeeze(mean(missLFP,1)),params);

% output data
output.hitLFP = hitLFP;
output.powerCWThit = powerCWThit;
output.tfhit = tfhit;
output.fwt = fwt;
output.frex = frex;
output.pnts = pnts;
output.Shit = Shit;

output.Sf = fhit;
output.Serrhit = Serrhit;
output.missLFP = missLFP;
output.powerCWTmiss = powerCWTmiss;
output.tfmiss = tfmiss;
output.Smiss = Smiss;
output.Serrmiss = Serrmiss;
end

function [trials,win] = makeSpikeWin(Spikes,spikeId,Fs)
%% Begin analysis across clusters
winLen = 4;
win = double(0:winLen:max(Spikes.SpikeSamples)/Fs);
trials = cell(1,length(win));
% win(1:400) = [];
for i = 1:length(win)-1
    spiketm = Spikes.Clusters(spikeId).spikeTime-win(i);
    idx = find(spiketm>=0 & spiketm < winLen);
    trials{i} = zeros(1,length(0:0.001:winLen));
    spiked = discretize(spiketm(idx),0:0.001:winLen);
    if ~isnan(spiked)
        trials{i}(spiked) = 1;
    end
end
end

function plotspikeCoherence(spikea,spikeb,spikeCoherence)

interval = 0:2000;
spiker = spikea;

figure;
%****************
disp('Plotting rastergrams (slow) ...');
subplot('position',[0.1 0.4 0.4 0.55]);
make_nice_spike_raster(spiker);
V = axis;
axis([0 size(spikea{1},2) V(3) V(4)]);
grid on;
ylabel('Trial Number');
title(fprintf('Unit rasters'));
%*************
subplot('position',[0.1 0.1 0.4 0.3]);
smooth_window = 25;  % give sigma of 12.5ms
make_nice_mean_raster(spiker,smooth_window);
V = axis;
axis([0 size(spikea{1},2) V(3) V(4)]);
plot([interval(1),interval(1)],[V(3),V(4)],'k-'); hold on;
plot([interval(end),interval(end)],[V(3),V(4)],'k-'); hold on;
ylabel('Firing Rate');
xlabel('Time (ms)');

spiker = spikeb;
%****************
disp('Plotting rasters2 per trial (slow) ...');
subplot('position',[0.57 0.4 0.4 0.55]);
make_nice_spike_raster(spiker);
V = axis;
axis([0 size(spikeb{1},2) V(3) V(4)]);
grid on;
ylabel('Trial Number');
title(fprintf('Unit rasters'));
%*************
subplot('position',[0.57 0.1 0.4 0.3]);
smooth_window = 25;  % give sigma of 12.5ms
make_nice_mean_raster(spiker,smooth_window);
V = axis;
axis([0 size(spikeb{1},2) V(3) V(4)]);
plot([interval(1),interval(1)],[V(3),V(4)],'k-'); hold on;
plot([interval(end),interval(end)],[V(3),V(4)],'k-'); hold on;
ylabel('Firing Rate');
xlabel('Time (ms)');

figure;
colo = 'rbgy';
CNUM=1;
subplot(2,2,1);
for ii = 1:CNUM
    H = semilogx(spikeCoherence.freq{ii},spikeCoherence.coho{ii},[colo(ii),'-']); hold on;
    set(H,'Linewidth',2);
    H = semilogx(spikeCoherence.freq{ii},(spikeCoherence.coho{ii}+spikeCoherence.scoho{ii}),[colo(ii),':']); hold on;
    H = semilogx(spikeCoherence.freq{ii},(spikeCoherence.coho{ii}-spikeCoherence.scoho{ii}),[colo(ii),':']);
    H = semilogx(spikeCoherence.freq{ii},spikeCoherence.rcoho{ii},[colo(ii),'--']); hold on;
    set(H,'Linewidth',1);
end
ylabel('Coherence Magnitude');
xlabel('Frequency (Hz)');
title(sprintf('Spike A with Spike B'));

subplot(2,2,2);
for ii = 1:CNUM
    H = plot(spikeCoherence.freq{ii},spikeCoherence.phaso{ii},[colo(ii),'-']); hold on;
    set(H,'Linewidth',2);
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.phaso{ii}+spikeCoherence.sphaso{ii}),[colo(ii),':']); hold on;
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.phaso{ii}-spikeCoherence.sphaso{ii}),[colo(ii),':']);
    %H = plot(freq{ii},rphaso{ii},[colo(ii),'--']); hold on;
    set(H,'Linewidth',1);
end
ylabel('Coherence Angle');
xlabel('Freq');
title(sprintf('Spike A with Spike B'));

subplot(2,2,3);
for ii = 1:CNUM
    H = plot(spikeCoherence.freq{ii},spikeCoherence.spika_pow{ii},[colo(ii),'-']); hold on;
    set(H,'Linewidth',2);
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.spika_pow{ii}+spikeCoherence.sspika_pow{ii}),[colo(ii),':']); hold on;
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.spika_pow{ii}-spikeCoherence.sspika_pow{ii}),[colo(ii),':']);
end
ylabel('Spike Unit A Pow');
xlabel('Frequency (Hz)');
title(sprintf('Spike A with Spike B'));

subplot(2,2,4);
for ii = 1:CNUM
    H = plot(spikeCoherence.freq{ii},spikeCoherence.spike_pow{ii},[colo(ii),'-']); hold on;
    set(H,'Linewidth',2);
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.spike_pow{ii}+spikeCoherence.sspike_pow{ii}),[colo(ii),':']); hold on;
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.spike_pow{ii}-spikeCoherence.sspike_pow{ii}),[colo(ii),':']);
end
ylabel('Spike Unit B Pow');
xlabel('Frequency (Hz)');
title(sprintf('Spike B with Spike A'));

end
%****************** function to make a nice raster plot ****************
function smorate = make_nice_mean_raster(spmat,smooth_window)
%*********** spmat1 and spmat2 are spike matrices of two conditions you wish to compare
%*********** smooth_window ... gaussian smoothing in millisecs
numconds = size(spmat,2);
if (numconds==2)
    colo = [[1,0,0];[0,0,1]];
else
    colo = [[1,0,0];[0,0,1];[0,1,0];[1,1,0]];
end

for k = 1:numconds
    spud = spmat{1,k};
    numtrials(k) = size(spud,1);
    smorate = gauss_smooth(sum( spud(1:numtrials(k),:))/....
        numtrials(k),smooth_window)*1000;
    [mSmo,idx] = max(smorate);
    if idx<1400 && idx>200
    smorate = gaussmf(1:length(smorate),[mSmo idx])'*mSmo;
    end
%     smorate = (smorate-min(smorate))/(max(smorate)-min(smorate));
%     smorate = zscore(smorate);
    H = plot(smorate,'k-'); hold on;
    set(H,'Color',colo(k,:));
end

end


%******************* make a rastergram of the actual spikes on each trial
function make_nice_spike_raster(spmat)

colo = [[1,1,1];[0,0,0]];
colormap(colo);

totspike = [];
for cubo = 1:size(spmat,2)
    totspike = [totspike; (cubo*spmat{1,cubo}) ];
    totspike = [totspike; (cubo*spmat{1,cubo}) ];
    
end
totspike = totspike + 1;
imagesc(totspike);
V = axis;
axis([V(1) V(2) 0 size(totspike,1)]);

end

%**************************************************************
function output = gauss_smooth(input, window)
% Smoothing function:
% output = smooth(input, window)
% "Window" is the total kernel width.
% Input array must be one-dimensional.

input_dims = ndims(input);
input_size = size(input);
if input_dims > 2 | min(input_size) > 1,
    disp('Input array is too large.');
    return
end

if input_size(2) > input_size(1),
    input = input';
    toggle_dims = 1;
else
    toggle_dims = 0;
end

if window/2 ~= round(window/2),
    window = window + 1;
end
halfwin = window/2;

input_length = length(input);
%********* gauss window +/- 1 sigma
x = -halfwin:1:halfwin;
kernel = exp(-x.^2/(window/2)^2);
kernel = kernel/sum(kernel);

padded(halfwin+1:input_length+halfwin) = input;
padded(1:halfwin) = ones(halfwin, 1)*input(1);
padded(length(padded)+1:length(padded)+halfwin) = ones(halfwin, 1)*input(input_length);

output = conv(padded, kernel);
output = output(window:input_length+window-1);

if toggle_dims == 1,
    output = output';
end
end