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
parameters.windowBeforePull = 3.5; % in seconds
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
IntanBehaviour.AvgHitTrace = mean(IntanBehaviour.AvgHitTrace,1);
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
Spikes = rejectSpikes(Spikes,0.4,0.25,IntanBehaviour.parameters); % Reject spikes here for further analysis
IntanBehaviour.reactionTime = 0;
[Spikes] = sortSpkLever(Spikes,IntanBehaviour);
[fpath,name,exts] = fileparts(ds_filename);
sessionName = [fpath,'/','Spikes.mat'];
% save(sessionName,"IntanBehaviour","fpath","parameters","-v7.3");
save(sessionName,"Spikes","IntanBehaviour","fpath","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
disp('Saved!')
%%
SqSpikes = zeros(size(Spikes.PSTH.hit.spks{1},1),size(Spikes.PSTH.hit.spks{1},2),length(Spikes.PSTH.hit.spks));
for n = 1:size(SqSpikes,3)
    SqSpikes(:,:,n) = Spikes.PSTH.hit.spks{n};
end
pullIndex = Spikes.PSTH.hit.pl;
pull1 = pullIndex(:,1);
pull2 = pullIndex(:,2);
pull3 = pullIndex(:,3);
tmin = 0;
tmax = size(SqSpikes,2)-1;

% Trial IDs, spiketimes, and neuron_ids are a flatten matrix of the 3d
% array. Trial IDs and spike times should have zero indexing
trial_ids = [];
spiketimes = [];
neuron_ids = [];
for ntrials = 1:size(SqSpikes,1)
    spkTemp = squeeze(SqSpikes(ntrials,:,:))';
    [r,c,v] = find(spkTemp==1);
    trial_ids = [trial_ids;(ntrials-1)*ones(size(r))];
    spiketimes = [spiketimes;c];
    neuron_ids = [neuron_ids;r];
end
neuron_ids = neuron_ids-1; % zero indexing
spiketimes = spiketimes-1; % zero indexing
[fpath,name,exts] = fileparts(ds_filename);
sessionName = [fpath,'\','spikes_to_Warp.mat'];
save(sessionName,"tmin","tmax","pull1","pull2","pull3","trial_ids","spiketimes","neuron_ids","fpath");
disp('Data Saved for warping')
%% Plot out spikes aligned to the pull response
% From the behaviour figure where we have the pull counts
%% Plot trial sorted by earliest first pull
% than by earliest second pull
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
    ylim([1 290])
    xline(0,'k','reward')
    spkTemp = squeeze(SqSpikes(:,:,neuron));
    spkTemp(spkTemp==0) = NaN;
    for n = 1:size(trialMask,1)
        scatter(time,n*spkTemp(firstPull(n),:),0.5,'filled','MarkerFaceColor',[0,0,0]),hold on
    end
    ylim([1 290])
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
[fpath,fname] = fileparts(ds_filename);
load(fullfile(fpath,'warpedSpks'))
% Plot out warped pulls
warpedSpks = getAlignedSqpulls(Spikes,warpedSpks,IntanBehaviour);
% Build spike rasters based on aligned data
%% Get statistics
% I want to plot out the peak z scored value at the time of pull for each
% unit. To do this, we should calculate calculate the residuals of each
% unit as a funciton of the pulls. That will tell us which neurons is
% tuned to the response... I think
% Parameters
% Parameters
bin_size = 20; % example bin size - adjust to your data/time resolution
analysis_window = [-0.5, 0.5]; % analysis window in binned time indices (adjust as needed)
all_pulls = warpedSpks.warpSpikes;
% Get sizes from original data
[num_trials, num_timebins, num_neurons, num_pulls] = size(all_pulls);

% Calculate the number of bins after binning
n_time_binned = floor(num_timebins / bin_size);
binTime = linspace(warpedSpks.warpedTime(1,3),warpedSpks.warpedTime(end,3),n_time_binned);
% Preallocate binned data array
binned_all_pulls = zeros(n_time_binned, num_neurons, num_pulls);

% Apply binning per neuron and pull
for p = 1:num_pulls
    for n = 1:num_neurons
        data = squeeze(all_pulls(:, :, n, p)); % trials × time
        binned_data = sum(getBin(data,bin_size))*(1000/bin_size);  % trials × binned_time
        binned_all_pulls(:, n, p) = binned_data;
    end
end

mean_baseline_neuron = squeeze(mean(binned_all_pulls(binTime<-2,:,:),[1,3]));
%%
% Preallocate mean response matrix: neurons × pulls
mean_responses = zeros(num_neurons, num_pulls);
wTime = warpedSpks.warpedTime;
binnedSpk = [];
% Calculate mean firing rate within analysis window per neuron/pull
% We want to normalize by the mean baseline
for pull = 1:num_pulls
    win = find(wTime(:,pull)>=analysis_window(1) & wTime(:,pull)<=analysis_window(2));
    data = squeeze(all_pulls(:, win, :, pull)); % trials × bins × neurons
    % Here we calculate the binned z score response of the neurons
    for neuron = 1:num_neurons
        spkTemp = squeeze(data(:,:,neuron));
        binnedSpk(:,:,neuron) = getBin(spkTemp,bin_size);
    end
    avg_window = sum(binnedSpk,1);    % average over time bin dimension -> trials × 1 × neurons
    avg_window = squeeze(avg_window)*(1000/bin_size); % trials × neurons

    mean_responses(:, pull) = mean(avg_window, 1); % mean over trials, result is 1 × neurons
end
%%
%here we subtract the baseline response to calculate the residuals so we
%now know how responsive the neuron was to the stimulus
residuals = mean_responses-mean_baseline_neuron'; 
modulationIndex = (mean_responses-mean_baseline_neuron')./(mean_responses+mean_baseline_neuron');
selectivity_index = zeros(num_neurons, 1);
best_pull = zeros(num_neurons, 1);

for n = 1:num_neurons
    responses = abs(residuals(n, :));
    [R_best, idx_best] = max(responses);
    R_other = mean(responses(setdiff(1:num_pulls, idx_best)));
    responseTot(n,:) = responses;
    selectivity_index(n) = (R_best - R_other) / (R_best + R_other);
    best_pull(n) = idx_best;
end

%% Plot Modulation for each pull

bin_edges = -1:0.1:1; % Adjust bin size if desired
colors = [0 0 1; 0 0.5 0; 0.7 0 0]; % blue, green, red for pulls
figure; hold on;

for p = 1:num_pulls
    subplot(1,3,p),
    % Bin and normalize
    counts = histcounts(modulationIndex(:,p), bin_edges, 'Normalization', 'probability');
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
temp = nan(length(selectivity_index),3);
y = [];
for p = 1:num_pulls
    sel_for_pull = selectivity_index(best_pull == p); % Example: your selectivity indices for this pull
    temp(1:length(sel_for_pull),p) = sel_for_pull;
    y(p) = length(sel_for_pull)/length(selectivity_index);
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


%%
function binned_data = getBin(data,bin_size)
[n_trials, n_time] = size(data);

% Trim time dimension to multiple of bin_size
n_time_trim = floor(n_time / bin_size) * bin_size;
data_trim = data(:, 1:n_time_trim);

% Reshape to [n_trials, bin_size, n_time_trim/bin_size]
data_reshaped = reshape(data_trim', bin_size, [], n_trials); 
% Note the transpose is so time is first dimension for reshaping

% Sum or average within bins (along first dimension)
binned_data = squeeze(mean(data_reshaped, 1))';  
% Output size: [n_trials, n_time_trim/bin_size]
end

