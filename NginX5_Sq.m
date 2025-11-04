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
%%
%%% Run Kilosort3 
% load only neccessary variables from memory mapped file
data = matfile(ds_filename);
fpath = data.fpath;
%Kilosort264FTestcode
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
parameters.effortPerturbation = 0
parameters.delay = 0.5; %reward delay
parameters.Fs = 1000; % Eventual downsampled data
parameters.ts = 1/parameters.Fs;
parameters.IntanFs = data.targetedFs;
parameters.rows = 64;
parameters.cols = 1;

[Behaviour] = readLeverSq(parameters,lfpTime,fname);
[IntanBehaviour] = readLeverIntan(parameters,data.amplifierTime,data.analogChannels(1,:),data.digitalChannels,Behaviour,1);

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
for i=1:IntanBehaviour.nCueHit
    plot(0:3000,smoothdata(IntanBehaviour.cueHitTrace(i).trace),'Color',[0 0 0 0.2],'LineWidth',1.5);
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
<<<<<<< Updated upstream
=======
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
    bin_centers = bin_edges(1:end-1) + diffs(bin_edges)/2;
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
    bin_centers = bin_edges(1:end-1) + diffs(bin_edges)/2;
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
Spikes.GPFA.HitEffort.dat = [Spikes.GPFA.MIHit.dat,Spikes.GPFA.effortperturb.dat];
for n = 1:length(IntanBehaviour.MIHitTrace)%+1:length(Spikes.GPFA.BaselineOpto.dat) %fix trials
    Spikes.GPFA.HitEffort.dat(n).trialId = n;
end
%%%
addpath(genpath('C:\Users\khan332\Documents\GitHub\NeuralTraj'));
addpath(genpath('mat_results'));
if exist('mat_results','dir'),rmdir('mat_results','s'),end
[Spikes.GPFA.resultHit,Spikes.GPFA.seqTrainHit] = gpfaAnalysis(Spikes.GPFA.hit.dat,1); %Run index
[Spikes.GPFA.resultMiss,Spikes.GPFA.seqTrainMiss] = gpfaAnalysis(Spikes.GPFA.miss.dat,2); %Run index
[Spikes.GPFA.resultMIHit,Spikes.GPFA.seqTrainMIHit] = gpfaAnalysis(Spikes.GPFA.MIHit.dat,3); %Run index
% [Spikes.GPFA.resultMIFA,Spikes.GPFA.seqTrainMIFA] = gpfaAnalysis(Spikes.GPFA.MIFA.dat,4); %Run index
[Spikes.GPFA.resultHitMiss,Spikes.GPFA.seqTrainHitMiss] = gpfaAnalysis(Spikes.GPFA.HitMiss.dat,5); %Run index
% [Spikes.GPFA.resultMIHitFA,Spikes.GPFA.seqTrainMIHitFA] = gpfaAnalysis(Spikes.GPFA.MIHitFA.dat,6); %Run index
[Spikes.GPFA.resultHitEffort,Spikes.GPFA.seqTrainHitEffort] = gpfaAnalysis(Spikes.GPFA.HitEffort.dat,7); %Run index
close all
sessionName = [fpath,'\','Spikes.mat'];
save(sessionName,"Spikes","IntanBehaviour","fpath","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
disp('Saved!')
%% For warping data only
[neuralDynamics] = getGPFASq(warpedSpks.pull1A.warpSpikes  ,IntanBehaviour);
%%
figure
plot(squeeze(neuralDynamics.hitOnly.X(2,:,:)),'color',[0.0 0.0 0.0 0.1]),hold on
plot(squeeze(mean(neuralDynamics.hitOnly.X(2,:,:),3)),'linewidth',2)
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
% Grab dimension of the data
x = squeeze(M1neuralDynamics.hiteffort.X(1,:,:));
y = squeeze(M1neuralDynamics.hiteffort.X(2,:,:));
z = squeeze(M1neuralDynamics.hiteffort.X(3,:,:));
% Now we can walk along time for all the trajectory points together
figure,hold on
for n = 1:size(x,2)
plot3(x(:,n),y(:,n),z(:,n),'color',[0 0 0 0.4])
plot3(x(1,n), y(1,n), z(1,n), 'o', 'MarkerFaceColor', [0.5 0.5 0.9], 'MarkerEdgeColor', 'k');
end
%%
x = squeeze(M1neuralDynamics.effort.X(1,:,:));
y = squeeze(M1neuralDynamics.effort.X(2,:,:));
z = squeeze(M1neuralDynamics.effort.X(3,:,:));
figure,hold on
for n = 1:80
plot3(x(:,n),y(:,n),z(:,n),'color',[0 0 0 0.4])
plot3(x(1,n), y(1,n), z(1,n), 'o', 'MarkerFaceColor', [0.9 0.5 0.5], 'MarkerEdgeColor', 'k');
end
%%
figure;
subplot(121),plot(squeeze(M1neuralDynamics.effort.X(2,:,:)))
subplot(122),plot(squeeze(M1neuralDynamics.hiteffort.X(2,:,:)))
%%
% Assuming x, y, z are [time x trials] matrices as per your code
x = horzcat(squeeze(M1neuralDynamics.hiteffort.X(1,:,:)),squeeze(M1neuralDynamics.effort.X(1,:,:)));
y = horzcat(squeeze(M1neuralDynamics.hiteffort.X(2,:,:)),squeeze(M1neuralDynamics.effort.X(2,:,:)));
z = horzcat(squeeze(M1neuralDynamics.hiteffort.X(3,:,:)),squeeze(M1neuralDynamics.effort.X(3,:,:)));
timeEnd = 250;
nTrials = size(x, 2);

v = VideoWriter('D:\SQLever\neural_trajectoriesM1Day15.avi'); % Name your output file
v.FrameRate = 20; % Set the frame rate
open(v);

figure('Color', 'w');
hold on
axis tight
view(30,45)
xlabel('X')
ylabel('Y')
zlabel('Z')
hitTrials = size(M1neuralDynamics.hiteffort.X,3);
for t = 1:timeEnd
    clf; % Clear the figure each frame
    hold on
    % Plot each trajectory up to time t

    plot3(x(1:t,:), y(1:t,:), z(1:t,:), 'Color', [0 0 0 0.4]);
    % Plot a dot for the current time point
    plot3(x(t,1:hitTrials), y(t,1:hitTrials), z(t,1:hitTrials), 'o', 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', 'k');
    plot3(x(t,hitTrials+1:end), y(t,hitTrials+1:end), z(t,hitTrials+1:end), 'o', 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', 'k');

    title(['Neural trajectories up to time = ', num2str(t)]);
    axis([min(x(:)), max(x(:)), min(y(:)), max(y(:)), min(z(:)), max(z(:))]);
    view(30, 45);
    grid on
    drawnow
    % Capture the frame and write to video
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);
%% Get time dependant initial conditions

[c,allTrials] = sort_hit_effort(IntanBehaviour);
x = horzcat(squeeze(M1neuralDynamics.hiteffort.X(1,:,:)),squeeze(M1neuralDynamics.effort.X(1,:,:)));
y = horzcat(squeeze(M1neuralDynamics.hiteffort.X(2,:,:)),squeeze(M1neuralDynamics.effort.X(2,:,:)));
z = horzcat(squeeze(M1neuralDynamics.hiteffort.X(3,:,:)),squeeze(M1neuralDynamics.effort.X(3,:,:)));
% sort trials
x = x(:,c);
y = y(:,c);
z = z(:,c);
timeEnd = 250;
nTrials = size(x, 2);

v = VideoWriter('D:\SQLever\neural_conditionsM1Day15.avi'); % Name your output file
v.FrameRate = 10; % Set the frame rate
open(v);
initial_azimuth = 30;
elevation = 45;

figure('Color', 'w');
hold on
axis tight
view(initial_azimuth,elevation)
xlabel('X')
ylabel('Y')
zlabel('Z')
hitTrials = size(M1neuralDynamics.hiteffort.X,3);
for t = 50
    clf; % Clear the figure each frame
    hold on
    % Plot each trajectory up to time t
    for trial = 1:size(x,2)
        plot3(x(t,1:trial), y(t,1:trial), z(t,1:trial), 'Color', [0 0 0 0.4]);
        % Plot a dot for the current time point based on effort
        if allTrials(2,trial)==1
            plot3(x(t,trial), y(t,trial), z(t,trial), 'o', 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', 'k');
        else
            plot3(x(t,trial), y(t,trial), z(t,trial), 'o', 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', 'k');
        end

        title(['Initial conditions on trial = ', num2str(trial)]);
        axis([min(x(:)), max(x(:)), min(y(:)), max(y(:)), min(z(:)), max(z(:))]);
        % Calculate current azimuth angle for rotation
        current_azimuth = mod(initial_azimuth + 50*trial, 360);
        view(current_azimuth, elevation);
        grid on
        drawnow
        % Capture the frame and write to video
        frame = getframe(gcf);
        writeVideo(v, frame);
    end
end

close(v);
%% Plot initial condition data
t = 50;
[c,allTrials] = sort_hit_effort(IntanBehaviour);
x = horzcat(squeeze(M1neuralDynamics.hiteffort.X(1,:,:)),squeeze(M1neuralDynamics.effort.X(1,:,:)));
y = horzcat(squeeze(M1neuralDynamics.hiteffort.X(2,:,:)),squeeze(M1neuralDynamics.effort.X(2,:,:)));
z = horzcat(squeeze(M1neuralDynamics.hiteffort.X(3,:,:)),squeeze(M1neuralDynamics.effort.X(3,:,:)));
figure,hold on
for trial = 1:size(x,2)
        % Plot a dot for the current time point based on effort
        if allTrials(2,trial)==1
            plot3(x(t,trial), y(t,trial), z(t,trial), 'o', 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', 'k');
        else
            plot3(x(t,trial), y(t,trial), z(t,trial), 'o', 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', 'k');
        end
end
view(30,30)
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
grid on
axis square
%% Perform K-means clustering
% Load your data (replace `initCond` with actual variable if different)
% initCond = ... % n x 3 matrix
[c,allTrials] = sort_hit_effort(IntanBehaviour);
x = horzcat(squeeze(M1neuralDynamics.hiteffort.X(1,:,:)),squeeze(M1neuralDynamics.effort.X(1,:,:)));
y = horzcat(squeeze(M1neuralDynamics.hiteffort.X(2,:,:)),squeeze(M1neuralDynamics.effort.X(2,:,:)));
z = horzcat(squeeze(M1neuralDynamics.hiteffort.X(3,:,:)),squeeze(M1neuralDynamics.effort.X(3,:,:)));
% sort trials
x = x(45,c);
y = y(45,c);
z = z(45,c);
initCond = [x',y',z'];
% Choose the number of clusters, e.g., 2 or 3 (can be tuned or assessed later)
numClusters = 3;
[idx, C] = kmeans(initCond, numClusters, 'Replicates', 10);
% Assume idx is your cluster assignment vector (values 1~3 for three clusters)
cluster_colors = [1 0.3 0.3;    % red
                  0.8 .3 0.8;    % blue
                  0.0 0.45 1];   % purple

figure;
hold on
for k = 1:3
    scatter3(initCond(idx==k,1), ...
             initCond(idx==k,2), ...
             initCond(idx==k,3), ...
             50, cluster_colors(k,:), 'filled');
end

% Overlay cluster centers (assuming C is centers array)
scatter3(C(:,1), C(:,2), C(:,3), 100, 'k', 'x', 'LineWidth', 2);
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
grid on
hold off
view(30,30)
axis square

%%% Calculate Mahalanobis distance of clusters
n = size(initCond, 1);
mahalDists = zeros(n, numClusters);
for k = 1:numClusters
    mu = C(k, :);
    sigma = cov(initCond(idx==k, :));      % Covariance of cluster k
    for i = 1:n
        mahalDists(i, k) = sqrt((initCond(i,:) - mu) / sigma * (initCond(i,:) - mu)'); % Mahalanobis distance formula
    end
end
% Visualize or output as needed
[sortedDist, trialOrder] = sort(mahalDists(:,1), 'descend');
mahalDists = mahalDists(trialOrder, :); % Now rows are ordered by their dist to Cluster 1
nTrials = size(mahalDists,1);


addpath(genpath('C:\Users\khan332\Documents\GitHub\slanCM'));
figure,hold on
plotNiceBars(mahalDists)
ylim([0 ceil(max(mahalDists,[],'all'))])
xlabel('Cluster');
ylabel('Mahalanobis Distance from Center');
title('Within-cluster Mahalanobis Distances');
normDist = (sortedDist - min(sortedDist)) / (max(sortedDist) - min(sortedDist));
cmap = slanCM('RdBu',nTrials); % or any other colormap
trialColors = cmap(round(normDist * (size(cmap,1)-1))+1, :);
colormap(slanCM('RdBu'))
c = colorbar;
c.Label.String = 'Sorted Mahalanobis Distance';
c.Label.FontSize = 8;
% Optionally set ticks to match the real value range:
c.Ticks = [0 0.5 1];
c.TickLabels = {num2str(min(sortedDist)), num2str(mean(sortedDist)), num2str(max(sortedDist))};
axis square
%% Silohette calculation
n = size(initCond,1);
pairwiseMahal = zeros(n,n);

% Use overall covariance for simplicity
Sigma = cov(initCond);

for i = 1:n
    for j = 1:n
        diffs = initCond(i,:) - initCond(j,:);
        pairwiseMahal(i,j) = sqrt(diffs / Sigma * diffs');
    end
end
silo = zeros(n,1);

for i = 1:n
    myCluster = idx(i);
    sameInds = find(idx == myCluster & (1:n)' ~= i);
    otherClusters = setdiff(unique(idx), myCluster);
    
    % a(i): mean Mahalanobis distance to same cluster
    a_i = mean(pairwiseMahal(i, sameInds));
    
    % b(i): minimum mean Mahalanobis distance to other clusters
    b_i = inf;
    for k = otherClusters'
        kInds = find(idx == k);
        b_ik = mean(pairwiseMahal(i, kInds));
        if b_ik < b_i
            b_i = b_ik;
        end
    end
    
    silo(i) = (b_i - a_i) / max(a_i, b_i);
end

% Plot it
figure;
h = histogram(silo, 11,'Normalization', 'probability', 'FaceColor', [0.5 0.8 1], 'EdgeColor','none');
set(gca, ...
    'TickDir', 'out', ...
    'Box', 'off', ...
    'FontSize', 14, ...
    'LineWidth', 1.5);
axis square;
xlabel('Silhouette Value', 'FontSize', 16);
ylabel('Probability', 'FontSize', 16);
title('Silhouette', 'FontSize', 16, 'FontWeight', 'normal');
% Set consistent limits for clarity
xlim([-0.5, 0.75]); % Adjust as needed for your data
ylim([0, max(h.Values)*1.1]);
% Compute and plot mean
m = median(silo);
skewness(silo)
yl = ylim;
hold on;
xline(m, '--k', ['Median = ' num2str(m, '%.2f')], ...
    'LineWidth', 2, ...
    'LabelOrientation', 'horizontal', ...
    'LabelHorizontalAlignment', 'center', ...
    'LabelVerticalAlignment', 'top', ...
    'FontSize', 8, ...
    'Color', [0.3 0.3 0.3]);
hold off;


s = silhouette(initCond,idx);
figure;
h = histogram(s, 'Normalization', 'probability', 'FaceColor', [0.5 0.8 1], 'EdgeColor','none');
set(gca, ...
    'TickDir', 'out', ...
    'Box', 'off', ...
    'FontSize', 14, ...
    'LineWidth', 1.5);
axis square;
xlabel('Silhouette Value', 'FontSize', 16);
ylabel('Probability', 'FontSize', 16);
title('Silhouette', 'FontSize', 16, 'FontWeight', 'normal');
% Set consistent limits for clarity
xlim([-0.2, 1]); % Adjust as needed for your data
ylim([0, max(h.Values)*1.1]);
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
time = linspace(-IntanBehaviour.parameters.windowBeforePull,IntanBehaviour.parameters.windowAfterPull,Spikes.GPFA.seqTrainHit(1).T);

% figure,hold on
% for n = 1:200
%     plot(time,squeeze(M1neuralDynamics.hit.speed.speed(1,:,n)),'color',[0 0 0 0.8])
%     for p = 1:length(IntanBehaviour.hitTrace(n).pullCount)
%         xline((IntanBehaviour.hitTrace(n).pullCount(p)-3500)/1000,'r')
%     end
% end

colors = [12,188,187;183,13,180]/255;
speedTotBaseline = smoothdata(squeeze(M1neuralDynamics.hiteffort.speed.speed(1,2:end,:)),1,'gaussian',10);
figure,hold on
plot(time(2:end),mean(speedTotBaseline,2),'color',colors(2,:),'linewidth',2),hold on
plot(time(2:end),mean(speedTotBaseline,2)+std(speedTotBaseline,[],2)/(sqrt(size(speedTotBaseline,2))),'color',colors(2,:),'linewidth',2)
plot(time(2:end),mean(speedTotBaseline,2)-std(speedTotBaseline,[],2)/(sqrt(size(speedTotBaseline,2))),'color',colors(2,:),'linewidth',2)
set(gca,'tickdir','out'),box off, axis square
xlim([-3.5,1.5])
%%
%%% LOCAL FUNCTIONS


%%
function plotNiceBars(totData)
means = nanmean(totData);          % Bar heights
sems = nanstd(totData) ./ sqrt(size(totData,1));   % Error bar (standard error)
b = bar(means, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k'); % Gray bars with black edge

% Overlay error bars
errorbar(1:size(totData,2), means, sems, 'k', 'LineStyle', 'none', 'LineWidth', 1);

% Overlay individual jittered points
xjitter = randn(size(totData))*0.01; % Controls point jitter
for i = 1:size(totData,2)
    scatter(i + xjitter(:,min(i,2)), totData(:,i), 18, 'o', ...
        'MarkerEdgeColor', [0.25 0.25 0.25], ...
        'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4);
end

nTrials = size(totData,1);
cm = slanCM('RdBu',nTrials); % Or use your favorite colormap
% Draw paired lines between columns 1 and 2
for j = 1:size(totData,1)
    if size(totData,2) >= 3  % If there are at least 3 columns
        xvals = [1 + xjitter(j,1), 2 + xjitter(j,2), 3 + xjitter(j,3)];
        yvals = [totData(j,1),    totData(j,2),    totData(j,3)];
        plot(xvals, yvals, '-', 'Color', cm(j,:), 'LineWidth', 1);
    else % Connect just columns 1 and 2
        xvals = [1 + xjitter(j,1), 2 + xjitter(j,2)];
        yvals = [totData(j,1),    totData(j,2)];
        plot(xvals, yvals, '-', 'Color', cm(j,:), 'LineWidth', 1);
    end
end

% Style similar to image
set(gca, 'XTick', 1:size(totData,2), 'XTickLabel', {'Second Pull', 'Third Pull', 'Polymer', 'Late'}, ...
    'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
ylabel('IPI (s)');
ylim([0 2]);

hold off;

%%% RUN STATS
[p, tbl, stats] = anova1(totData, [], 'off'); % columns as groups
results = multcompare(stats, 'Display', 'off'); % Pairwise comparisons

disp(['ANOVA p-value: ', num2str(p)]);
alpha = 0.05; % significance level
sigPairs = results(results(:,6) < alpha, :); % rows where p < 0.05
hold on;
ylims = ylim;

% vertical height offset for significance lines above bars
baseY = max(means + sems) * 1.05;  
offsetStep = max(means + sems) * 0.05; 
if all(results(:,6) >= 0.05) % No significant pairwise differences
    % Extract F statistic from ANOVA table
    Fstat = cell2mat(tbl(2,5)); % Assumes standard anova1 output tbl
    p_anova = p;
    % Place text on plot upper corner
    xPos = size(totData,2)/2;
    yPos = max(means + sems) * 2.4;
    text(xPos, yPos, sprintf('ANOVA F=%.2f, p=%.3f', Fstat, p_anova), ...
        'HorizontalAlignment', 'left', 'FontSize', 10);
    % Add pairwise stars or p-values as before (your existing code)
end

for i = 1:size(sigPairs,1)
    x1 = sigPairs(i,1);
    x2 = sigPairs(i,2);
    y = baseY + (i-1)*offsetStep;
    
    % Draw line connecting bars
    plot([x1 x1 x2 x2], [y y+offsetStep y+offsetStep y], 'k-', 'LineWidth', 1);
    
    % Add star above the line
    text(mean([x1 x2]), y + offsetStep*0.1, '*', 'HorizontalAlignment', 'center', ...
        'FontSize', 16, 'FontWeight', 'bold');
end
hold off;
end
>>>>>>> Stashed changes

