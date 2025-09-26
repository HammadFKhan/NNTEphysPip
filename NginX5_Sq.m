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
%% Run Kilosort3 
% load only neccessary variables from memory mapped file
data = matfile(ds_filename,'Writable',true);
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
parameters.perturbEffort = 0; %perturb effort after first pull

parameters.windowBeforePull = 3.5; % in seconds % technically window before reward
parameters.windowAfterPull = 1.5; % in seconds
parameters.windowBeforeCue = 1.5; % in seconds
parameters.windowAfterCue = 3.5; % in seconds
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
% temp1 = arrayfun(@(x) x.LFPtime(1), IntanBehaviour.cueHitTrace);
% temp1 = vertcat(temp1,ones(1,IntanBehaviour.nCueHit)); %  write 1 for reward given
% temp2 = arrayfun(@(x) x.LFPtime(1), IntanBehaviour.cueMissTrace);
% temp2 = vertcat(temp2,zeros(1,IntanBehaviour.nCueMiss)); %  write 0 for no reward given
% temp = [temp1,temp2];
% [~,idx] = sort(temp(1,:)); %sort by occurance
% IntanBehaviour.ITI = temp(:,idx);
IntanBehaviour.parameters = parameters;
Behaviour.parameters = parameters;
[Behaviour.cleanedpullCounts, Behaviour.pullIndices] = cleanTimeoutSequences(Behaviour,0);

%% Plot behaviour
figure
for i=1:length(IntanBehaviour.hitTrace)
    plot(0:5000,smoothdata(IntanBehaviour.MIHitTrace(i).trace),'Color',[0 0 0 0.2],'LineWidth',1.5);
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
IntanBehaviour.AvgHitTrace = mean(IntanBehaviour.AvgHitTrace,1);
figure,
for n = 1:length(IntanBehaviour.missTrace)
    plot(0:5000,smoothdata(IntanBehaviour.missTrace(n).trace),'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
    missTrace(n,:) = IntanBehaviour.missTrace(n).trace;
end
IntanBehaviour.AvgMissTrace = mean(missTrace,1);
IntanBehaviour.AvgHitTrace = mean(IntanBehaviour.AvgHitTrace,1);
%%% sequence dynamics
pullIndex = vertcat(IntanBehaviour.hitTrace.pullCount);
figure,imagesc(Behaviour.cleanedpullCounts)
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
addpath(genpath('C:\Users\khan332\Documents\GitHub\Kilosort')) % path to kilosort folder
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
Spikes = leverPSTHSq(Spikes,IntanBehaviour);
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
Spikes = rejectSpikes(Spikes,0.1,0.15,IntanBehaviour.parameters); % Reject spikes here for further analysis
IntanBehaviour.reactionTime = 0;
[Spikes] = sortSpkLever(Spikes,IntanBehaviour);
[fpath,name,exts] = fileparts(ds_filename);
sessionName = [fpath,'/','Spikes.mat'];
% save(sessionName,"IntanBehaviour","fpath","parameters","-v7.3");
save(sessionName,"Spikes","IntanBehaviour","fpath","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
disp('Saved!')
%% Prep data for warping
prepWrap(Spikes,ds_filename)
%% Plot out spikes aligned to the pull response
% From the behaviour figure where we have the pull counts
%% Plot trial sorted by earliest first pull
% than by earliest second pull
pullIndex = vertcat(IntanBehaviour.hitTrace.pullCount);
[ft,firstPull] = sort(pullIndex(:,1));
[sc,secondPull] = sort(pullIndex(:,2));
[trialMask] = getAUTOResponse(IntanBehaviour,1); %expFlag = 1
figure,
time = linspace(-IntanBehaviour.parameters.windowBeforePull,IntanBehaviour.parameters.windowAfterPull,size(trialMask,2));
figure
color = [46,49,149]/255;
for n = 1:size(trialMask,1)
    scatter(time,n*trialMask(n,:),10,'filled','MarkerFaceColor',color),hold on
end
xlim([-2.5, 0.5])
ylim([1 290])
xline(0,'k','reward')
axis square

time = linspace(-IntanBehaviour.parameters.windowBeforePull,IntanBehaviour.parameters.windowAfterPull,size(trialMask,2));
figure,subplot(121)
color = [46,49,149;]/255;
for n = 1:size(trialMask,1)
    scatter(time,n*trialMask(firstPull(n),:),10,'filled','MarkerFaceColor',color),hold on
end
xlim([-2.5, 0.5])
ylim([1 290])
xline(0,'k','reward')
axis square
subplot(122),
color = [46,149,49]/255;
for n = 1:size(trialMask,1)
    scatter(time,n*trialMask(secondPull(n),:),10,'filled','MarkerFaceColor',color),hold on
end
xlim([-2.5, 0.5])
ylim([1 290])
xline(0,'k','reward')
axis square
%% sort to first pull
% Nice neurons 1 3 6 7 10 19 17
time = linspace(-IntanBehaviour.parameters.windowBeforePull,IntanBehaviour.parameters.windowAfterPull,size(trialMask,2));
figure
count = 1;
color = [46,49 149]/255;
for neuron = [3 4 23 24]
    subplot(3,2,count)
    xlim([-2.5, 0.5])
    ylim([1 100])
    xline(0,'k','reward')
    spkTemp = squeeze(warpedSpks.pull3A.warpSpikes(:,:,neuron));
    spkTemp(spkTemp==0) = NaN;
    for n = 1:size(trialMask,1)
        scatter(time,n*spkTemp(firstPull(n),:),0.5,'filled','MarkerFaceColor',[0,0,0]),hold on
    end
    ylim([1 100])
    xline(0,'k','reward')
    xlim([-2.5, 0.5])
    title(['Neuron ' num2str(neuron)])
    count = count+1;
    for n = 1:size(trialMask,1)
        scatter(time(ft(n)),n*trialMask(firstPull(n),ft(n)),5,'filled','MarkerEdgeColor','none','MarkerFaceColor',color),hold on
    end
end
%% sort to second pull
% Nice neurons 1 3 6 7 10 19 17
time = linspace(-IntanBehaviour.parameters.windowBeforePull,IntanBehaviour.parameters.windowAfterPull,size(trialMask,2));
figure
count = 1;
color = [46,149,49]/255;
for neuron = 1:30
    subplot(6,5,count)
    xlim([-2.5, 0.5])
    ylim([1 100])
    xline(0,'k','reward')
    spkTemp = Spikes.PSTH.hit.spks{neuron};
    spkTemp(spkTemp==0) = NaN;
    for n = 1:size(trialMask,1)
        scatter(time,n*spkTemp(secondPull(n),:),5,'filled','MarkerFaceColor',[0,0,0]),hold on
    end
    ylim([1 100])
    xline(0,'k','reward')
    xlim([-2.5, 0.5])
    for n = 1:size(trialMask,1)
        scatter(time(sc(n)),n*trialMask(secondPull(n),sc(n)),10,'filled','MarkerEdgeColor',color,'MarkerFaceColor','none'),hold on
    end
    count = count+1;
end
%% Load in warped data and make psth based on warping models
% 
if ~exist('fpath','var')
[fpath,fname] = fileparts(ds_filename);
end
load(fullfile(fpath,'warpedSpks'))
% Plot out warped pulls
warpedSpks = getAlignedSqpulls(Spikes,warpedSpks,IntanBehaviour);
close all
% Build spike rasters based on aligned data
sessionName = [fpath,'/','warpedSpks.mat'];
% save(sessionName,"IntanBehaviour","fpath","parameters","-v7.3");
save(sessionName,"warpedSpks","Spikes","IntanBehaviour","fpath","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
disp('Saved!')
%% Get statistics
warpedSpks = getWarpSpkStats(warpedSpks);
%% Plot Modulation for each pull
bin_edges = -1:0.1:1; % Adjust bin size if desired
colors = [0 0 1; 0 0.5 0; 0.7 0 0]; % blue, green, red for pulls
figure; hold on;
num_pulls = size(warpedSpks.stats.modulationIndex,2);
for p = 1:num_pulls
    subplot(1,3,p),
    % Bin and normalize
    counts = histcounts(warpedSpks.stats.modulationIndex(:,p), bin_edges, 'Normalization', 'probability');
    bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
    % Stairs plot
    stairs(bin_centers, counts, 'Color', colors(p,:), 'LineWidth', 2);
    xlabel('Modulation Index');
    ylabel('Probability');
    title('Modulation Index');
    hold off;
    set(gca,'tickdir','out'),axis square
end


%% temp(1:length(sel_for_pull),p) = sel_for_pull;
temp = nan(length(warpedSpks.stats.selectivityIndex),3);
y = [];
for p = 1:num_pulls
    sel_for_pull = warpedSpks.stats.selectivityIndex(warpedSpks.stats.bestPull == p); % Example: your selectivity indices for this pull
    temp(1:length(sel_for_pull),p) = sel_for_pull;
    y(p) = length(sel_for_pull)/length(warpedSpks.stats.selectivityIndex);
end

color = [46,49,179;46,149,49;179,49,46]/255;
figure,violinplot(temp,[],'ViolinColor',color);
xlabel('Pull Preference (Selectivity) Index')
ylabel('Number of Neurons')
legend({'Pull 1', 'Pull 2', 'Pull 3'}, 'Location', 'Best')
title('Distribution of Neuronal Pull Preference')
hold off;box off,set(gca,'tickdir','out'),axis square
figure,
x = [1];
bar(x,y,'stacked')
hold off;box off,set(gca,'tickdir','out'),axis square
xlim([0 2])
%%
% 
figure;
subplot(131),scatter(responseTot(:,1), responseTot(:,3), 24, 'filled','k');
hold on;
plot([0 1], [0 1], 'k--', 'LineWidth', 1.2); % y = x reference line
hold off;
xlabel('Mean Spike Response', 'FontWeight', 'bold');
ylabel('Mean Spike Response', 'FontWeight', 'bold');
title('Pull 1 vs Pull 3', 'FontWeight', 'bold');
set(gca, 'FontSize', 8, 'Box', 'off', 'GridAlpha', 0.4);
axis square;
box off,set(gca,'tickdir','out')

subplot(132),scatter(responseTot(:,1), responseTot(:,2), 24, 'filled','k');
hold on;
plot([0 1], [0 1], 'k--', 'LineWidth', 1.2); % y = x reference line
hold off;
xlabel('Mean Spike Response', 'FontWeight', 'bold');
ylabel('Mean Spike Response', 'FontWeight', 'bold');
title('Pull 1 vs Pull 2', 'FontWeight', 'bold');
set(gca, 'FontSize', 8, 'Box', 'off', 'GridAlpha', 0.4);
axis square;
box off,set(gca,'tickdir','out')
subplot(133),scatter(responseTot(:,2), responseTot(:,3), 24, 'filled','k');
hold on;
plot([0 1], [0 1], 'k--', 'LineWidth', 1.2); % y = x reference line
hold off;
xlabel('Mean Spike Response', 'FontWeight', 'bold');
ylabel('Mean Spike Response', 'FontWeight', 'bold');
title('Pull 2 vs Pull 3', 'FontWeight', 'bold');
set(gca, 'FontSize', 8, 'Box', 'off', 'GridAlpha', 0.4);
axis square;
box off,set(gca,'tickdir','out')
%%
bin_edges = 0:0.05:1; % Adjust bin size if desired
colors = [0 0 1; 0 0.5 0; 0.7 0 0]; % blue, green, red for pulls
temp = nan(length(selectivity_index),3);
figure; hold on;

for p = 1:num_pulls
    sel_for_pull = selectivity_index(best_pull == p); % Example: your selectivity indices for this pull
    temp(1:length(sel_for_pull),p) = sel_for_pull;
    % Bin and normalize
    counts = histcounts(sel_for_pull, bin_edges, 'Normalization', 'count');
    bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
    % Stairs plot
    stairs(bin_centers, counts, 'Color', colors(p,:), 'LineWidth', 2);
end

xlabel('Pull Preference (Selectivity) Index');
ylabel('Probability');
legend({'Pull 1', 'Pull 2', 'Pull 3'}, 'Location', 'Best');
title('Normalized Pull Preference Index Distribution (Stairs Plot)');
hold off;
%% Neural Trajectory Segementation using GPFA
% Note that we concatenate trial conditions as to apply the same models for
% statistical comparison (ie. hit vs miss, hit vs FA, opto vs no opto)
Spikes = makeSpikeGPFA(Spikes);
Spikes.GPFA.HitMiss.dat = [Spikes.GPFA.hit.dat,Spikes.GPFA.miss.dat];
for n = 1:length(IntanBehaviour.hitTrace)%+1:length(Spikes.GPFA.BaselineOpto.dat) %fix trials
    Spikes.GPFA.HitMiss.dat(n).trialId = n;
end
Spikes.GPFA.MIHitFA.dat = [Spikes.GPFA.MIHit.dat,Spikes.GPFA.MIFA.dat];
for n = length(IntanBehaviour.MIHitTrace)+1:length(Spikes.GPFA.MIHitFA.dat) %fix trials
    Spikes.GPFA.MIHitFA.dat(n).trialId = n;
end
%%%
addpath(genpath('C:\Users\khan332\Documents\GitHub\NeuralTraj'));
addpath(genpath('mat_results'));
if exist('mat_results','dir'),rmdir('mat_results','s'),end
[Spikes.GPFA.resultHit,Spikes.GPFA.seqTrainHit] = gpfaAnalysis(Spikes.GPFA.hit.dat,1); %Run index
[Spikes.GPFA.resultMiss,Spikes.GPFA.seqTrainMiss] = gpfaAnalysis(Spikes.GPFA.miss.dat,2); %Run index
[Spikes.GPFA.resultMIHit,Spikes.GPFA.seqTrainMIHit] = gpfaAnalysis(Spikes.GPFA.MIHit.dat,3); %Run index
[Spikes.GPFA.resultMIFA,Spikes.GPFA.seqTrainMIFA] = gpfaAnalysis(Spikes.GPFA.MIFA.dat,4); %Run index
[Spikes.GPFA.resultHitMiss,Spikes.GPFA.seqTrainHitMiss] = gpfaAnalysis(Spikes.GPFA.HitMiss.dat,5); %Run index
[Spikes.GPFA.resultMIHitFA,Spikes.GPFA.seqTrainMIHitFA] = gpfaAnalysis(Spikes.GPFA.MIHitFA.dat,6); %Run index
close all
sessionName = [fpath,'\','Spikes.mat'];
save(sessionName,"Spikes","IntanBehaviour","fpath","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
disp('Saved!')

%% Neural Trajectory Analysis
%IntanBehaviour.parameters = parameters;
%neuralTrajAnalysis(Spikes,Waves1,IntanBehaviour1);
[M1neuralDynamics,M1waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);
%%
time = linspace(-IntanBehaviour.parameters.windowBeforePull,IntanBehaviour.parameters.windowAfterPull,Spikes.GPFA.seqTrainHit(1).T);
figure,
for n = 1:length(Spikes.GPFA.seqTrainHit)
    subplot(311),plot(time,Spikes.GPFA.seqTrainHit(n).xorth(1,:),'Color',[0 0 0 0.4]),hold on,box off,set(gca,'fontsize',18)
    subplot(312),plot(time,Spikes.GPFA.seqTrainHit(n).xorth(2,:),'Color',[0 0 0 0.4]),hold on,box off,set(gca,'fontsize',18)
    subplot(313),plot(time,Spikes.GPFA.seqTrainHit(n).xorth(3,:),'Color',[0 0 0 0.4]),hold on,box off,set(gca,'fontsize',18)
end
%%
pullIndex = vertcat(IntanBehaviour.hitTrace.pullCount);
x = squeeze(M1neuralDynamics.hit.X(1,:,:));
y = squeeze(M1neuralDynamics.hit.X(2,:,:));
z = squeeze(M1neuralDynamics.hit.X(3,:,:));
figure,hold on
for n = 1:286
plot3(x(:,n),y(:,n),z(:,n),'color',[0 0 0 0.4])
end
%%
colors = [12,188,187;183,13,180]/255;
timeIndex = linspace(1,size(x,1),5001);

xp = mean(x,2);
yp = mean(y,2);
zp = mean(z,2);
figure
plot3(xp,yp,zp,'color',colors(2,:)),hold on
pIm = floor(mean(pullIndex));
pIm = [1 pIm 3499]; %add reward 
for n = 1:length(pIm)
    id = floor(timeIndex(pIm(n)));
    scatter3(xp(id),yp(id),zp(id),20,'filled'),hold on
end

%%
% figure,hold on
% for n = 1:200
%     plot(time,squeeze(M1neuralDynamics.hit.speed.speed(1,:,n)),'color',[0 0 0 0.8])
% %     for p = 1:length(IntanBehaviour.hitTrace(n).pullCount)
% %         xline((IntanBehaviour.hitTrace(n).pullCount(p)-3500)/1000,'r')
% %     end
% end
colors = [12,188,187;183,13,180]/255;
speedTotBaseline = smoothdata(squeeze(M1neuralDynamics.hit.speed.speed(1,2:end,:)),1,'gaussian',10);
figure,hold on
plot(time(2:end),mean(speedTotBaseline,2),'color',colors(2,:),'linewidth',2),hold on
plot(time(2:end),mean(speedTotBaseline,2)+std(speedTotBaseline,[],2)/(sqrt(size(speedTotBaseline,2))),'color',colors(2,:),'linewidth',2)
plot(time(2:end),mean(speedTotBaseline,2)-std(speedTotBaseline,[],2)/(sqrt(size(speedTotBaseline,2))),'color',colors(2,:),'linewidth',2)
set(gca,'tickdir','out'),box off, axis square
xlim([-3.5,1.5])
