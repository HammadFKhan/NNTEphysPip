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
parameters.opto = 1; % 1 - opto ON , 0 - opto OFF
parameters.cool = 0; % No Cool 
parameters.windowBeforePull = 3; % in seconds
parameters.windowAfterPull = 2; % in seconds
parameters.windowBeforeCue = 1.5; % in seconds
parameters.windowAfterCue = 1.5; % in seconds
parameters.windowBeforeMI = 1.5; % in seconds 
parameters.windowAfterMI = 3.5; % in seconds 
parameters.perturbEffort = 0;
parameters.delay = 0.5; %reward delay
parameters.Fs = 1000; % Eventual downsampled data
parameters.ts = 1/parameters.Fs;
parameters.IntanFs = data.targetedFs;
parameters.rows = 64;
parameters.cols = 1;
lfpTime = data.amplifierTime;
[Behaviour] = readLeverSq(parameters,lfpTime);
[IntanBehaviour] = readLeverIntanSq(parameters,data.amplifierTime,data.analogChannels(1,:),data.digitalChannels,Behaviour,0);

IntanBehaviour.reactionTime = arrayfun(@(x) x.pullCount(3)-x.pullCount(1),IntanBehaviour.hitTrace)/1000;
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
Spikes = rejectSpikes(Spikes,0.25,0.25,IntanBehaviour.parameters); % Reject spikes here for further analysis
[Spikes] = sortSpkLever(Spikes,IntanBehaviour);
[fpath,name,exts] = fileparts(ds_filename);
sessionName = [fpath,'/','Spikes.mat'];
% save(sessionName,"IntanBehaviour","fpath","parameters","-v7.3");
save(sessionName,"Spikes","IntanBehaviour","fpath","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
disp('Saved!')
%% Prep data for warping
prepWrap(Spikes,ds_filename)%% Behaviour
%%
optoTrials = vertcat(IntanBehaviour.hitTrace.optoFlag);   % 0/1 labels
rtHit       = abs(vertcat(IntanBehaviour.reactionTime));
noOptoHitrt = rtHit(optoTrials == 0);
optoHitrt   = rtHit(optoTrials == 1);

m_rt  = [mean(noOptoHitrt,'omitnan'), mean(optoHitrt,'omitnan')];
sem_rt = [std(noOptoHitrt,'omitnan')/sqrt(numel(noOptoHitrt)), ...
          std(optoHitrt,'omitnan')  /sqrt(numel(optoHitrt))];

col_noOpto = [0.4 0.4 0.4];        % gray base
col_pts_no = [0.2 0.2 0.2];        % dark gray points
col_pts_op = [0 123 167]/255;      % cerulean points [web:115][web:116]

figure; hold on

% Bars (semi‑transparent gray)
b = bar(1:2, m_rt, 'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
b.FaceAlpha = 0.5;

% Error bars
errorbar(1:2, m_rt, sem_rt, 'k', 'LineStyle','none', 'LineWidth',1);

% Overlay data points with horizontal jitter
jitter = 0.08;

x1 = 1 + (rand(size(noOptoHitrt))-0.5)*2*jitter;
x2 = 2 + (rand(size(optoHitrt))-0.5)*2*jitter;

scatter(x1, noOptoHitrt, 25, col_pts_no, 'filled', 'MarkerFaceAlpha',0.8);
scatter(x2, optoHitrt,   25, col_pts_op, 'filled', 'MarkerFaceAlpha',0.8);

set(gca,'XTick',1:2,'XTickLabel',{'No opto','Opto'});
ylabel('Reaction time (s)');
set(gca,'Box','off','TickDir','out','FontSize',12);


% ---------- stats ----------
% Non‑parametric between-groups test (rank‑sum / Mann‑Whitney U) 
[p,~,stats] = ranksum(noOptoHitrt, optoHitrt);   

fprintf('Rank-sum: p = %.3g, z = %.3f\n', p, stats.zval);

% ---------- plot (assuming you already drew bars + points) ----------
hold on
yl = ylim;
txt = sprintf('p = %.3g', p);
text(1.5, yl(2)*0.95, txt, 'HorizontalAlignment','center', ...
     'FontWeight','bold');   
%% Plot out behaviour
% Extract lever traces: [time x nTrials]
hitTrace = horzcat(IntanBehaviour.MIHitTrace.trace);    % concatenate cells
optoTrials = vertcat(IntanBehaviour.hitTrace.optoFlag);   % 0/1 labels

% Split trials
trace_noOpto = hitTrace(:, optoTrials == 0);   % [time x nNo]
trace_opto   = hitTrace(:, optoTrials == 1);   % [time x nOp]

t = (1:size(hitTrace,1));   % or actual time vector in ms/s

% Mean ± SEM for each condition
mean_noOpto = mean(trace_noOpto, 2, 'omitnan');
mean_opto   = mean(trace_opto,   2, 'omitnan');

sem_noOpto = std(trace_noOpto, 0, 2, 'omitnan') ./ sqrt(size(trace_noOpto,2));
sem_opto   = std(trace_opto,   0, 2, 'omitnan') ./ sqrt(size(trace_opto,2));

% Trial‑wise difference (opto − no opto) for overlapping trials count
nPairs = min(size(trace_noOpto,2), size(trace_opto,2));
diff_traces = trace_opto(:,1:nPairs) - trace_noOpto(:,1:nPairs);  % [time x nPairs]
mean_diff = mean(diff_traces, 2, 'omitnan');
sem_diff  = std(diff_traces, 0, 2, 'omitnan') ./ sqrt(nPairs);

figure;

% -------- Top: mean lever traces --------
subplot(2,1,1); hold on
col_no = [0.4 0.4 0.4];
col_op = [0 123 167]/255;   % cerulean

% shaded SEM (simple patch)
fill([t fliplr(t)], [(mean_noOpto-sem_noOpto)' fliplr((mean_noOpto+sem_noOpto)')], ...
     col_no, 'FaceAlpha',0.2, 'EdgeColor','none');
fill([t fliplr(t)], [(mean_opto-sem_opto)'   fliplr((mean_opto+sem_opto)')], ...
     col_op, 'FaceAlpha',0.2, 'EdgeColor','none');

plot(t, mean_noOpto, 'Color', col_no, 'LineWidth',1.5);
plot(t, mean_opto,   'Color', col_op, 'LineWidth',1.5);

ylabel('Lever position');
title('Lever traces: No opto vs Opto');
set(gca,'Box','off','TickDir','out');
xlim([0 3000])
% -------- Bottom: difference trace (opto − no opto) --------
subplot(2,1,2); hold on
fill([t fliplr(t)], [(mean_diff-sem_diff)' fliplr((mean_diff+sem_diff)')], ...
     col_op, 'FaceAlpha',0.2, 'EdgeColor','none');
plot(t, mean_diff, 'Color', col_op, 'LineWidth',1.5);

yline(0,'k--');
xlabel('Time');
ylabel('\Delta lever (opto - no opto)');
set(gca,'Box','off','TickDir','out');
xlim([0 3000])
%% Plot out example tagged units
neuronId  = 2;% 6 9 12 16
tempSpk   = Spikes.PSTH.MIHit.spks{neuronId};   % [nTrials x nTime]
optoTrials = vertcat(IntanBehaviour.hitTrace.optoFlag);   % [nTrials x 1], 0/1

nonOptoIdx = optoTrials == 0;
optoIdx = find(optoTrials == 1);

spk_nonOpto = tempSpk(nonOptoIdx,:);   % non‑opto trials

tempSpk   = Spikes.PSTH.MIHit.spks{neuronId};      % [nTrials x nTime], 0/1

spk_mod = tempSpk;        % copy to modify
optoAd      = 0;          % 0 = off, 1 = add opto-driven spikes
if optoAd==1
    spk_opto = adOpto(tempSpk,optoIdx);
else
    spk_opto    = tempSpk(optoIdx,:);      % opto trials
end
t = (1:size(tempSpk,2));               % time axis (samples or ms)

figure;

% -------- 1) Raster: non-opto --------
subplot(2,2,1); hold on
[row,col] = find(spk_nonOpto);
scatter(t(col), row, 6, [0.5 0.5 0.5], 'filled');                 % [web:6]
ylabel('Trials (no opto)');
title(['Example neuron: ', num2str(neuronId)]);

set(gca,'YDir','reverse','Box','off','TickDir','out');
ylim([0 20])
% -------- 2) Raster: opto --------
subplot(2,2,2); hold on
[row,col] = find(spk_opto);
scatter(t(col), row, 6, [0 123 167] / 255, 'filled');         % magenta for opto
ylabel('Trials (opto)');
title('Optogenetic trials');
set(gca,'YDir','reverse','Box','off','TickDir','out');

% -------- 3) Mean rate: non-opto --------
binSize = 20;                           % samples or ms per bin
edges = 1:binSize:(size(tempSpk,2)+1);
centers = edges(1:end-1) + binSize/2;

% non-opto PSTH
cnt = zeros(size(spk_nonOpto,1), numel(edges)-1);
for b = 1:numel(edges)-1
    idx = edges(b):edges(b+1)-1;
    cnt(:,b) = sum(spk_nonOpto(:,idx),2);
end
rate_nonOpto = mean(cnt,1) * (1000/binSize);            % Hz [web:19]

subplot(2,2,3); hold on
plot(centers, rate_nonOpto, 'Color',[0.5 0.5 0.5],'LineWidth',1.5);
xlabel('Time (ms)');
ylabel('Rate (spks/s)');
title('Mean rate (no opto)');
set(gca,'Box','off','TickDir','out');

% -------- 4) Mean rate: opto --------
cnt = zeros(size(spk_opto,1), numel(edges)-1);
for b = 1:numel(edges)-1
    idx = edges(b):edges(b+1)-1;
    cnt(:,b) = sum(spk_opto(:,idx),2);
end
rate_opto = mean(cnt,1) * (1000/binSize);               % Hz [web:19]

subplot(2,2,4); hold on
plot(centers, rate_opto, 'Color',[0 123 167] / 255,'LineWidth',1.5);
xlabel('Time (ms)');
ylabel('Rate (spks/s)');
title('Mean rate (opto)');
set(gca,'Box','off','TickDir','out');
% Add vertical lines at each pulse
optoStart = 1500;      % ms
pulseStep = 50;        % ms between pulses
tEnd      = 2500;      % end of plotting window
pulseTimes = optoStart:pulseStep:tEnd;
for pt = pulseTimes
    xline(pt, '-', 'Color', [0.6 0.6 0.6]);   % dotted gray lines [web:68][web:88]
end

xlabel('Time (ms)');
ylabel('Rate (spks/s)');
title('Mean rate (opto)');
set(gca,'Box','off','TickDir','out');

%% Calculate proportion of neurons that are opto tagged
nBoot   = 1000;
alpha   = 0.05;          % 95% CI
binSize = 20;
stim_ms = 1500;
stim_bin = ceil(stim_ms / binSize);

nNeurons = numel(Spikes.PSTH.MIHit.spks);
MI_all   = nan(1,nNeurons);
tagged   = false(1,nNeurons);

for neuronId = 1:nNeurons
    tempSpk = Spikes.PSTH.hit.spks{neuronId};
    if isempty(tempSpk), continue; end

    optoTrials  = vertcat(IntanBehaviour.hitTrace.optoFlag);
    spk_nonOpto = tempSpk(optoTrials==0,:);
    spk_opto    = tempSpk(optoTrials==1,:);

    % mean rates post‑stim (as before)
    edges = 1:binSize:(size(tempSpk,2)+1);
    nBins = numel(edges)-1;
    rB = getRate(spk_nonOpto,binSize,stim_bin);
    rO = getRate(spk_opto,binSize,stim_bin);
    den = rO + rB;
    if abs(den) < 1e-6
        MI_all(neuronId) = NaN;
        continue
    end
    MI_all(neuronId) = (rO - rB) / den;

    % ---- bootstrap CI on MI ----
    mi_boot = nan(1,nBoot);
    nB = size(spk_nonOpto,1);
    nO = size(spk_opto,1);
    for k = 1:nBoot
        idxB = randi(nB,[nB 1]);
        idxO = randi(nO,[nO 1]);
        rB_k = getRate(spk_nonOpto(idxB,:),binSize,stim_bin);
        rO_k = getRate(spk_opto(idxO,:),binSize,stim_bin);
        denk = rO_k + rB_k;
        if abs(denk) < 1e-6
            mi_boot(k) = NaN;
        else
            mi_boot(k) = (rO_k - rB_k) / denk;
        end
    end
    mi_boot = sort(mi_boot(~isnan(mi_boot)));
    if isempty(mi_boot), continue; end

    lo = mi_boot(round((alpha/2)*numel(mi_boot)));
    hi = mi_boot(round((1-alpha/2)*numel(mi_boot)));

    % neuron is "tagged" if CI does not include 0
    tagged(neuronId) = (lo > 0);
    minMI = 0.05;  % example threshold
    tagged(neuronId) = (lo > 0) && (MI_all(neuronId) > minMI);
    % old version
%     tagged(neuronId) = (lo > 0) || (hi < 0);
    fprintf('  Neuron %d: MI = %.3f, 95%% CI [%.3f, %.3f], tagged = %d\n', ...
        neuronId, MI_all(neuronId), lo, hi, tagged(neuronId));
end
Spikes.BiPOLES.tagged = tagged;
Spikes.BiPOLES.MI = MI_all;
fprintf('Finished.\n');
fprintf('Total tagged neurons: %d / %d (%.1f%%)\n', ...
        sum(tagged), nNeurons, 100*sum(tagged)/nNeurons);
sessionName = [fpath,'/','Spikes.mat'];
% save(sessionName,"IntanBehaviour","fpath","parameters","-v7.3");
save(sessionName,"Spikes","IntanBehaviour","fpath","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
disp('Saved!')

%% Plot proportions
m_rt  = sum(Spikes.BiPOLES.tagged)/length(Spikes.BiPOLES.tagged);


col_noOpto = [0.4 0.4 0.4];        % gray base
col_pts_no = [0.2 0.2 0.2];        % dark gray points
col_pts_op = [0 123 167]/255;      % cerulean points [web:115][web:116]

figure; hold on

% Bars (semi‑transparent gray)
b = bar(1, m_rt, 'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
b.FaceAlpha = 0.5;

% Error bars
% errorbar(1:2, m_rt, sem_rt, 'k', 'LineStyle','none', 'LineWidth',1);

% Overlay data points with horizontal jitter
% jitter = 0.08;
% 
% x1 = 1 + (rand(size(noOptoHitrt))-0.5)*2*jitter;
% x2 = 2 + (rand(size(optoHitrt))-0.5)*2*jitter;

% scatter(x1, noOptoHitrt, 25, col_pts_no, 'filled', 'MarkerFaceAlpha',0.8);
% scatter(x2, optoHitrt,   25, col_pts_op, 'filled', 'MarkerFaceAlpha',0.8);

set(gca,'XTick',1,'XTickLabel',{'Opto'});
ylabel('Tagged Neurons (%)');
set(gca,'Box','off','TickDir','out','FontSize',12);
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
%% Analyze neural dynamics of opto and non opto trials
% Note that we concatenate trial conditions as to apply the same models for
% statistical comparison (ie. hit vs miss, hit vs FA, opto vs no opto)
Spikes = makeSpikeGPFA(Spikes);
Spikes.GPFA.HitMiss.dat = [Spikes.GPFA.hit.dat,Spikes.GPFA.miss.dat];
for n = 1:length(Spikes.GPFA.hit.dat)+1:length(Spikes.GPFA.miss.dat) %fix trials
    Spikes.GPFA.HitMiss.dat(n).trialId = n;
end
Spikes.GPFA.MIHitFA.dat = [Spikes.GPFA.MIHit.dat,Spikes.GPFA.MIFA.dat];
for n = length(Spikes.GPFA.MIHit.dat)+1:length(Spikes.GPFA.MIHitFA.dat) %fix trials
    Spikes.GPFA.MIHitFA.dat(n).trialId = n;
end
%%%
addpath(genpath('C:\Users\khan332\Documents\GitHub\NeuralTraj'));
addpath(genpath('mat_results'));
if exist('mat_results','dir'),rmdir('mat_results','s'),end
SpikesNonOpto = Spikes.GPFA.hit.dat;
optoTrials = vertcat(IntanBehaviour.hitTrace.optoFlag);   % [nTrials x 1]
nonOptoIdx = optoTrials == 0;
optoIdx    = optoTrials == 1;
SpikesNonOpto = SpikesNonOpto(nonOptoIdx);
[Spikes.GPFA.resultHit,Spikes.GPFA.seqTrainHit] = gpfaAnalysis(Spikes.GPFA.hit.dat,1); %Run index
% the hit only has the non opto, we can check to see how much it changes
% the loadings
[Spikes.GPFA.resultHitnoOpto,Spikes.GPFA.seqTrainHitnoOpto] = gpfaAnalysis(SpikesNonOpto,2); %Run index
[Spikes.GPFA.resultMiss,Spikes.GPFA.seqTrainMiss] = gpfaAnalysis(Spikes.GPFA.miss.dat,3); %Run index
[Spikes.GPFA.resultMIHit,Spikes.GPFA.seqTrainMIHit] = gpfaAnalysis(Spikes.GPFA.MIHit.dat,4); %Run index
[Spikes.GPFA.resultMIFA,Spikes.GPFA.seqTrainMIFA] = gpfaAnalysis(Spikes.GPFA.MIFA.dat,5); %Run index
[Spikes.GPFA.resultHitMiss,Spikes.GPFA.seqTrainHitMiss] = gpfaAnalysis(Spikes.GPFA.HitMiss.dat,6); %Run index
[Spikes.GPFA.resultMIHitFA,Spikes.GPFA.seqTrainMIHitFA] = gpfaAnalysis(Spikes.GPFA.MIHitFA.dat,7); %Run index
close all
sessionName = [fpath,'\','Spikes.mat'];
save(sessionName,"Spikes","IntanBehaviour","fpath","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
disp('Saved!')
%% Neural Trajectory Analysis
%IntanBehaviour.parameters = parameters;
%neuralTrajAnalysis(Spikes,Waves1,IntanBehaviour1);
[neuralDynamics,waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);

%%
optoTrials = vertcat(IntanBehaviour.hitTrace.optoFlag);   % [nTrials x 1]

nonOptoIdx = optoTrials == 0;
optoIdx    = optoTrials == 1;
pullIdx = [];
for n = 1:length(IntanBehaviour.MIHitTrace)
    pullIdx(n,:) = [IntanBehaviour.MIHitTrace(n).pullCount(1:3), IntanBehaviour.MIHitTrace(n).pullCount(3)+500];
end
pullIdxm = mean(pullIdx,1);
pullIdxm = floor(pullIdxm/20);
% nonOptoIdx = 1:40;
% optoIdx = 41:76;
% dim x time x trials  ->  time x trials per dim
x = squeeze(neuralDynamics.hitOnly.X(1,:,:));   % [time x trials]
y = squeeze(neuralDynamics.hitOnly.X(2,:,:));
z = squeeze(neuralDynamics.hitOnly.X(3,:,:));

xno = x(:,nonOptoIdx);   yno = y(:,nonOptoIdx);   zno = z(:,nonOptoIdx);
xo  = x(:,optoIdx);      yo  = y(:,optoIdx);      zo  = z(:,optoIdx);

% mean trajectories
mx_no = mean(xno,2)*2; my_no = mean(yno,2)*2; mz_no = mean(zno,2)*2;
mx_o  = mean(xo,2);  my_o  = mean(yo,2);  mz_o  = mean(zo,2);

% colors
col_no = [0.6 0.6 0.6];           % light gray baseline
col_o  = [0 123 167]/255;         % cerulean for opto [web:115]

reactionTime = vertcat(IntanBehaviour.reactionTime);
startIdx = 1;
stimIdx  = 75;                    % stimulus bin
rtIdx = floor((1500+mean(reactionTime(nonOptoIdx))*1000)/20);
figure; hold on
factorScale = (length(optoTrials)-sum(optoTrials))/sum(optoTrials);
% baseline trajectory (thin gray)
plot3(mx_no, my_no, mz_no, 'Color', col_no, 'LineWidth', 2);

% opto trajectory (thicker cerulean)
plot3(mx_o, my_o, mz_o, 'Color', col_o, 'LineWidth', 2.5);

% markers for start and stim non opto trajectory
plot3(mx_no(startIdx), my_no(startIdx), mz_no(startIdx), ...
      'o', 'MarkerSize', 8, 'MarkerFaceColor', col_no, 'MarkerEdgeColor','none');
for n = 1:length(pullIdxm)
    plot3(mx_no(pullIdxm(n)), my_no(pullIdxm(n)), mz_no(pullIdxm(n)), ...
      'o', 'MarkerSize', 8, 'MarkerFaceColor', col_no, 'MarkerEdgeColor','none');
end

plot3(mx_no(stimIdx),  my_no(stimIdx),  mz_no(stimIdx),  ...
      'o', 'MarkerSize', 8, 'MarkerFaceColor', col_no, 'MarkerEdgeColor','none');


plot3(mx_no(rtIdx),  my_no(rtIdx),  mz_no(rtIdx),  ...
      'o', 'MarkerSize', 8, 'MarkerFaceColor', col_no, 'MarkerEdgeColor','none');

% markers at start and stim along opto trajectory
plot3(mx_o(startIdx), my_o(startIdx), mz_o(startIdx), ...
      'o', 'MarkerSize', 8, 'MarkerFaceColor', col_o, 'MarkerEdgeColor','none');

for n = 1:length(pullIdxm)
plot3(mx_o(pullIdxm(n)),  my_o(pullIdxm(n)),  mz_o(pullIdxm(n)),  ...
      'o', 'MarkerSize', 8, 'MarkerFaceColor', col_o, 'MarkerEdgeColor','none');
end

% rtIdx = floor((1500+mean(reactionTime(optoIdx))*1000)/20);
% plot3(mx_o(rtIdx),  my_o(rtIdx),  mz_o(rtIdx),  ...
%       'o', 'MarkerSize', 8, 'MarkerFaceColor', col_o, 'MarkerEdgeColor','none');

set(gca,'Box','off','TickDir','out','XColor','k','YColor','k','ZColor','k');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
view(3);
%% Speed
optoTrials = vertcat(IntanBehaviour.hitTrace.optoFlag);
nonOptoIdx = optoTrials == 0;
optoIdx    = optoTrials == 1;

% neuralDynamics.hitOnly.speed.speed: dim x time x trials
spd = neuralDynamics.MIhit.speed.speed;   % assume dim 1 = speed

t = linspace(-1.5,3.5,size(spd,2));   % or construct time vector in s

spd_no = squeeze(spd(1,:,nonOptoIdx));   % [time x nNo]
spd_no = zscore(spd_no);
spd_o  = squeeze(spd(1,:,optoIdx));      % [time x nOp]
spd_o = zscore(spd_o);

mn_no  = mean(spd_no,2);
mn_o   = mean(spd_o,2);
sem_no = std(spd_no,0,2)./sqrt(size(spd_no,2));
sem_o  = std(spd_o,0,2)./sqrt(size(spd_o,2));

col_no = [0.4 0.4 0.4];
col_o  = [0 123 167]/255;   % cerulean [web:115][web:116]

cueTime   = 0;      % s
miTime1   = floor((mean(reactionTime(nonOptoIdx))*1000))/1000;   % first MI boundary
miTime2   = floor((mean(reactionTime(optoIdx))*1000))/1000;    % second MI boundary

figure; hold on

% shaded SEM: non‑opto
fill([t fliplr(t)], [(mn_no-sem_no)' fliplr((mn_no+sem_no)')], ...
     col_no, 'FaceAlpha',0.2, 'EdgeColor','none');
% shaded SEM: opto
fill([t fliplr(t)], [(mn_o-sem_o)' fliplr((mn_o+sem_o)')], ...
     col_o, 'FaceAlpha',0.2, 'EdgeColor','none');

% mean traces
plot(t, mn_no, 'Color', col_no, 'LineWidth',2);
plot(t, mn_o,  'Color', col_o,  'LineWidth',2);

% vertical lines
xline(cueTime, '--', 'Color',[0 0 0],   'LineWidth',1.5);   % cue
xline(miTime1,'--', 'Color',[0.5 0.5 0.5],'LineWidth',1.5); % MI window 1
xline(miTime2,'--', 'Color',col_o,      'LineWidth',1.5);   % MI window 2

xlabel('Time from cue (s)');
ylabel('Trajectory speed');
set(gca,'Box','off','TickDir','out','FontName','Helvetica','FontSize',10);
%% Plot PC loadings
% hit only
PCloadingswOpto = Spikes.GPFA.resultHit.kern.estParams.L;
PCloadingsnoOpto = Spikes.GPFA.resultHitnoOpto.kern.estParams.L;
figure,subplot(121),bar(PCloadingswOpto(:,1))
subplot(122),bar(PCloadingsnoOpto(:,1))
%% Cosine of loading vectors
% Row-normalize so each neuron vector has unit norm
E1 = PCloadingswOpto;
E2 = PCloadingsnoOpto;
[cosSim_neuron, meanCosSim] = neuronEmbeddingCosineSimilarity(E1, E2);

fprintf('Mean neuron-wise cosine similarity = %.3f\n', meanCosSim);

% For visualization:
figure; histogram(cosSim_neuron, 'BinWidth', 0.05);
xlabel('Cosine similarity of neuron embeddings');
ylabel('Count');
title('Similarity of neuron loadings across subspaces');
%% Now make a barplot based on loading of tagged and untagged neurons
TaggedNCosine = cosSim_neuron(Spikes.BiPOLES.tagged==1);
unTaggedNCosine = cosSim_neuron(Spikes.BiPOLES.tagged==0);
cosineDat = nan(max([length(TaggedNCosine),length(unTaggedNCosine)]),2);
cosineDat(1:length(TaggedNCosine),1) = TaggedNCosine;
cosineDat(1:length(unTaggedNCosine),2) = unTaggedNCosine;
figure,hold on
plotNiceBars(cosineDat)
%% Plot explained variance of tagged vs untagged neurons
% loadings matrix
[~, nPCs] = size(PCloadingsnoOpto);

% Build groups as index vectors
groups = cell(1,2);
groups{1} = find(Spikes.BiPOLES.tagged == 1);      % tagged
groups{2} = find(Spikes.BiPOLES.tagged == 0);     % untagged

groupNames = {'Tagged', 'Untagged'};

nGroups = numel(groups);

% Denominator: sum of squared loadings per PC (should be ~1)
colNormSq = sum(PCloadingsnoOpto.^2, 1);     % 1 x nPCs

% Fraction of squared loading mass per group, per PC
fracLoad = nan(nGroups, nPCs);

for g = 1:nGroups
    idx = groups{g};                             % neuron indices for this group
    num = sum(PCloadingsnoOpto(idx,:).^2, 1);    % 1 x nPCs
    fracLoad(g,:) = num ./ colNormSq;           % fraction per PC
end

% Choose PCs to show (e.g., first 6)
pcsToPlot = 1:min(6, nPCs);
dataToPlot = fracLoad(:, pcsToPlot)';   % nPCsShown x nGroups

% Colors (tagged vs untagged)
cBlue = [0 0.45 0.74];      % tagged
cGrey = [0.5 0.5 0.5];      % untagged
colors = [cBlue; cGrey];

figure; hold on;

hBar = bar(pcsToPlot, dataToPlot, 'grouped');
for g = 1:numel(hBar)
    set(hBar(g), ...
        'FaceColor', colors(g,:), ...
        'EdgeColor', 'none');
end

% Add white edges only on the outer outline if you like:
for g = 1:numel(hBar)
    hBar(g).LineWidth = 0.5;
end

set(gca, 'Box', 'off', ...
         'TickDir', 'out', ...
         'LineWidth', 1, ...
         'FontName', 'Helvetica', ...
         'FontSize', 10);

xlim([pcsToPlot(1)-0.5, pcsToPlot(end)+0.5]);
ylim([0 1]);

xlabel('Principal component');
ylabel('Fraction of squared loading');

xticks(pcsToPlot);
xticklabels(arrayfun(@(k) sprintf('%d', k), pcsToPlot, 'UniformOutput', false));

legend(groupNames, 'Location', 'northoutside', ...
       'Orientation', 'horizontal', ...
       'Box', 'off');

set(gca,'tickdir','out');
%% Project loadings into 3d space
% Choose which PCs to visualize
pcX = 1;   % x-axis PC
pcY = 2;   % y-axis PC
pcZ = 3;
Lx = PCloadingsnoOpto(:, pcX);
Ly = PCloadingsnoOpto(:, pcY);
Lz = PCloadingsnoOpto(:, pcZ);
% Logical masks
tagged   = find(Spikes.BiPOLES.tagged == 1);
untagged = find(Spikes.BiPOLES.tagged == 0);

% Colors
cBlue = [0 0.45 0.74];     % tagged
cGrey = [0.6 0.6 0.6];     % untagged

figure; hold on;

% Untagged neurons (background)
scatter3(Lx(untagged), Ly(untagged),Lz(untagged), 20, ...
        'MarkerFaceColor', cGrey, ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.5);

% Tagged neurons (highlight)
scatter3(Lx(tagged), Ly(tagged),Lz(tagged), 30, ...
        'MarkerFaceColor', cBlue, ...
        'MarkerEdgeColor', 'w', ...
        'LineWidth', 0.5);

axis equal;

set(gca, 'Box', 'off', ...
         'TickDir', 'out', ...
         'LineWidth', 1, ...
         'FontName', 'Helvetica', ...
         'FontSize', 10);

xlabel(sprintf('PC %d loading', pcX));
ylabel(sprintf('PC %d loading', pcY));

legend({'Untagged','Tagged'}, ...
       'Location', 'northoutside', ...
       'Orientation', 'horizontal', ...
       'Box', 'off');
%%
X = neuralDynamics.hitOnly.X;
% X: nPC x nTime x nTrials
[nPC, nTime, nTrials] = size(X);

t = linspace(-3,2,250);  % 1 x nTime, e.g. aligned to behavior (you should already have this)

% Mean across trials
PCmean = mean(X, 3);   % nPC x nTime

% Select PCs
pcsToPlot = 1:min(6, nPC);

cGrey = [0.5 0.5 0.5];
cBlue = [0 0.45 0.74];

figure;

for k = 1:numel(pcsToPlot)
    pcIdx = pcsToPlot(k);
    
    subplot(numel(pcsToPlot), 1, k); hold on;
    
    % Example: plot mean PC time course in blue
    plot(t, PCmean(pcIdx,:), 'Color', cBlue, 'LineWidth', 1.5);
    
    % Optional: add SEM shading if you want
    % PCsem = std(X(pcIdx,:,:), 0, 3) / sqrt(nTrials);
    % fill_between could be emulated with patch
    
    set(gca, 'Box', 'off', ...
             'TickDir', 'out', ...
             'LineWidth', 1, ...
             'FontName', 'Helvetica', ...
             'FontSize', 10);
    
    ylabel(sprintf('PC %d', pcIdx));
    
    if k == numel(pcsToPlot)
        xlabel('Time (ms)');
    else
        set(gca, 'XTickLabel', []);  % no x labels for upper panels
    end
    
    % Optionally: vertical line at event (e.g., stim onset)
    % xline(0, '--', 'Color', [0.3 0.3 0.3]);
end
%%
%%% LOCAL FUNCTIONS

function [cosSim_neuron, meanCosSim] = neuronEmbeddingCosineSimilarity(E1, E2)
% E1, E2: [nNeurons x nDims] neuron embedding matrices
% Returns:
%   cosSim_neuron: [nNeurons x 1] cosine similarity per neuron
%   meanCosSim: scalar, average across neurons

    % Basic checks
    if ~isequal(size(E1), size(E2))
        error('E1 and E2 must have the same size [nNeurons x nDims].');
    end

    % Flatten to double
    E1 = double(E1);
    E2 = double(E2);

    % Compute norms per neuron (row-wise)
    n1 = sqrt(sum(E1.^2, 2));   % [nNeurons x 1]
    n2 = sqrt(sum(E2.^2, 2));   % [nNeurons x 1]

    % Avoid division by zero: set zero-norm rows to eps
    n1(n1 == 0) = eps;
    n2(n2 == 0) = eps;

    % Dot product per neuron across dimensions
    dotProd = sum(E1 .* E2, 2);  % [nNeurons x 1]

    % Cosine similarity per neuron: cos(theta_i)
    cosSim_neuron = dotProd ./ (n1 .* n2);  % in [-1, 1]

    % Optional summary: mean cosine similarity across neurons
    meanCosSim = mean(cosSim_neuron, 'omitnan');
end

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
%%
function spk_opto = adOpto(tempSpk,optoIdx)
spk_mod = tempSpk;        % copy to modify


optoFreq    = 20;         % Hz
binSize_ms  = 1;          % your current spike bin size
optoStart   = 1500;       % first bin to consider (ms)
nPulses     = 20;          % how many pulses to simulate
stepBins    = round((1000/optoFreq)/binSize_ms);   % 50 bins

% Example: probability that can vary across pulses (bins)
% length must be >= nPulses
p_vec = linspace(0.9, 0.6, nPulses);   % low → high probability


for k = 1:numel(optoIdx)
    tr = optoIdx(k);

    pulse = 0;
    for t = optoStart:stepBins:size(spk_mod,2)
        pulse = pulse + 1;
        if pulse > numel(p_vec)
            break
        end

        p_this = p_vec(pulse);          % probability for this bin

        % draw spike for this bin in this trial
        if rand < p_this
            spk_mod(tr,t) = 1;
        end
    end
end


spk_opto    = spk_mod(optoIdx,:);      % opto trials
end

% helper to get post‑stim mean rate
function r = getRate(spk,binSize,stim_bin)
edges = 1:binSize:(size(spk,2)+1);
nBins = numel(edges)-1;
cnt = zeros(size(spk,1),nBins);
for b = 1:nBins
    idx = edges(b):edges(b+1)-1;
    cnt(:,b) = sum(spk(:,idx),2);
end
rate = mean(cnt,1) * (1000/binSize);
r = mean(rate(stim_bin:end));
end
