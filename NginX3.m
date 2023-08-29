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

%% LFP
set(0,'DefaultFigureWindowStyle','normal')
LFP = fastpreprocess_filtering(lfp,Fs);
LFP = bestLFP(LFP);
LFP = bandFilter(LFP,'depth'); % Extract LFPs based on 'depth' or 'single'
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
for i=1:Behaviour.nHit
    plot(0:2000,smoothdata(Behaviour.hitTrace(i).trace),'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
    try
        hitTrace(i,:) = smoothdata(Behaviour.hitTrace(i).rawtrace);
    catch
        continue
    end
end

%% CSD and spectrogram
leverLFPAnalysis(LFP,Behaviour)
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
function leverLFPAnalysis(LFP,Behaviour)
CSDoutputhit = [];waveletHit = [];waveletMiss = [];powerCWThit = [];CSDoutputmiss = []; hitLFP = [];missLFP = [];
eIdx = find(s.sorted_probe_wiring(:,2)==0);
linearProbe = LFP.medianLFP(eIdx,:);
params.tapers = [5 9];
params.Fs = 1000;
params.fpass = [0 80];
params.err = [2 0.05];
for i = 1:100
    hitWin = [Behaviour.hit(i,3)-1000, Behaviour.hit(i,3)+1000]; %intan position
    timestamphit(i,:) = [Behaviour.hit(i,4)-1, Behaviour.hit(i,4)+1]; %s
    hitLFP(:,:,i) = linearProbe(1:4,hitWin(1):hitWin(2));
    [waveletHit(:,:,i),fwavelet] = cwt(mean(hitLFP(:,:,i),1),1000,'FrequencyLimit',[4 80]);
    [powerCWThit(:,:,i), fwt] = calCWTSpectogram(mean(hitLFP(:,:,i),1),0:2000,1000,10,[4 80],0);
    [CSDoutputhit(:,:,i)]  = CSD(hitLFP(:,:,i)'/1E6,1000,20E-6);
end
for i = 1:Behaviour.nMiss
    missWin = [Behaviour.miss(i,3)-1000, Behaviour.miss(i,3)+1000]; %ms
    timestampmiss(i,:) = [Behaviour.miss(i,4)-1, Behaviour.miss(i,4)+1]; %s
    missLFP(:,:,i) = linearProbe(1:4,missWin(1):missWin(2));
    [waveletMiss(:,:,i),fwavelet] = cwt(mean(missLFP(:,:,i),1),1000,'FrequencyLimit',[4 80]);
    [powerCWTmiss(:,:,i), fwt] = calCWTSpectogram(mean(missLFP(:,:,i),1),0:2000,1000,10,[4 80],0);
    [CSDoutputmiss(:,:,i)]  = CSD(missLFP(:,:,i)'/1E6,1000,20E-6);
end
for i = 1:size(linearProbe,1)
    tfmiss(:,:,i) = itpc(linearProbe(i,:),timestampmiss,1000,1);
    tfhit(:,:,i) = itpc(linearProbe(i,:),timestamphit,1000,1);
end
figure,plot(smoothdata(squeeze(mean(tfhit,1)),'gaussian',50),'k')
hold on,plot(smoothdata(squeeze(mean(tfmiss,1)),'gaussian',50),'r')
[Shits,fhits,Serrhits]=mtspectrumc(squeeze(mean(hitLFP,1)),params);
[Smiss,fmiss,Serrmiss] = mtspectrumc(squeeze(mean(missLFP,1)),params);

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