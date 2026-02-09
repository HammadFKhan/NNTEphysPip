%% Bipoles spike analysis
%% Parameters for behaviour
data = matfile(ds_filename); % ds_filename comes from loadme.mat
% check if data directory matches where the file originated; if not we note
% the new directory path
parameters.experiment = 'cue'; % self - internally generated, cue - cue initiated
parameters.opto = 1; % 1 - opto ON , 0 - opto OFF
parameters.cool = 0; % No Cool 
parameters.windowBeforePull = 1.5; % in seconds
parameters.windowAfterPull = 1.5; % in seconds
parameters.windowBeforeCue = 1.5; % in seconds
parameters.windowAfterCue = 1.5; % in seconds
parameters.windowBeforeMI = 1.5; % in seconds 
parameters.windowAfterMI = 1.5; % in seconds 
parameters.Fs = 1000; % Eventual downsampled data
parameters.ts = 1/parameters.Fs;
parameters.IntanFs = data.targetedFs;
parameters.rows = 64;
parameters.cols = 1;
parameters.BiPOLES = 1;
if ~exist('fname','var')
    [enfile,enpath] = uigetfile('Y:\Hammad\Ephys\SeqProject\*.csv');
    if isequal(enfile,0)
        disp('User selected Cancel');
    else
        disp(['User selected ', fullfile(enpath,enfile)]);
    end
else
    [enpath,enfile,ext] = fileparts(fname);
    disp(['User selected ', fullfile(enpath,enfile)]);
    enfile = [enfile,ext];
end
[Behaviour] = readLeverBiPOLESv2(enpath,enfile,parameters,data.amplifierTime,0);
[IntanBehaviour] = readLeverIntanBiPOLESv2(parameters,data.amplifierTime,data.analogChannels(1,:),data.digitalChannels,Behaviour,1);
% Calculate ITI time for trials and reward/no reward sequence
temp1 = arrayfun(@(x) x.LFPtime(1), IntanBehaviour.cueHitTrace);
temp1 = vertcat(temp1,ones(1,IntanBehaviour.nCueHit)); %  write 1 for reward given
temp2 = arrayfun(@(x) x.LFPtime(1), IntanBehaviour.cueMissTrace);
temp2 = vertcat(temp2,zeros(1,IntanBehaviour.nCueMiss)); %  write 0 for no reward given
temp = [temp1,temp2];
[~,idx] = sort(temp(1,:)); %sort by occurance
IntanBehaviour.ITI = temp(:,idx);
%% Combine opto cue hit and baseline cue hit
IntanBehaviour.optoCueHitTrace = rmfield(IntanBehaviour.optoCueHitTrace, {'MIIndex', 'cueIndex'});
IntanBehaviour.cueHitTrace = rmfield(IntanBehaviour.cueHitTrace, 'rewardIndex');
A = IntanBehaviour.cueHitTrace;       % destination struct
B = IntanBehaviour.optoCueHitTrace;   % source struct

TA = struct2table(A);
TB = struct2table(B);

Tmerged = [TA; TB];        % or vertcat(TA, TB)

IntanBehaviour.cueHitTrace = table2struct(Tmerged);
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
%load chanMap64Sharp
load chanMap64M
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
%% Behaviour
optoTrials = vertcat(IntanBehaviour.cueHitTrace.opto);   % 0/1 labels
rtHit       = abs(vertcat(IntanBehaviour.cueHitTrace.reactionTime));
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
hitTrace = horzcat(IntanBehaviour.cueHitTrace.trace);    % concatenate cells
optoTrials = vertcat(IntanBehaviour.cueHitTrace.opto);   % 0/1 labels

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
neuronId  = 12;% 12 16
tempSpk   = Spikes.PSTH.hit.spks{neuronId};   % [nTrials x nTime]
optoTrials = vertcat(IntanBehaviour.cueHitTrace.opto);   % [nTrials x 1], 0/1

nonOptoIdx = optoTrials == 0;
optoIdx = find(optoTrials == 1);

spk_nonOpto = tempSpk(nonOptoIdx,:);   % non‑opto trials

tempSpk   = Spikes.PSTH.hit.spks{neuronId};      % [nTrials x nTime], 0/1
optoTrials = vertcat(IntanBehaviour.cueHitTrace.opto);   % 0/1

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
tEnd      = 3000;      % end of plotting window
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

nNeurons = numel(Spikes.PSTH.hit.spks);
MI_all   = nan(1,nNeurons);
tagged   = false(1,nNeurons);

for neuronId = 1:nNeurons
    tempSpk = Spikes.PSTH.hit.spks{neuronId};
    if isempty(tempSpk), continue; end

    optoTrials  = vertcat(IntanBehaviour.cueHitTrace.opto);
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
%% Refine tag neurons (optional)

allSpk      = Spikes.PSTH.hit.spks;                % 1 x nNeurons cell
optoTrials  = vertcat(IntanBehaviour.cueHitTrace.opto);   % 0/1
nonOptoIdx  = optoTrials == 0;
optoIdx     = find(optoTrials == 1);

nNeurons = numel(allSpk);
spk_nonOpto_all = cell(1,nNeurons);
spk_opto_all    = cell(1,nNeurons);

for neuronId = 1:nNeurons
    tempSpk = allSpk{neuronId};        % [nTrials x nTime]
    if isempty(tempSpk), continue; end

    % non‑opto trials are always original spikes
    spk_nonOpto_all{neuronId} = tempSpk(nonOptoIdx,:);

    if tagged(neuronId)
        % tagged: add opto‑locked spikes only in opto trials
        spk_opto_all{neuronId} = adOpto(tempSpk, optoIdx);
    else
        % untagged: keep original spikes in opto trials
        spk_opto_all{neuronId} = tempSpk(optoIdx,:);
    end
end

allSpk     = Spikes.PSTH.hit.spks;                 % 1 x nNeurons cell
optoTrials = vertcat(IntanBehaviour.cueHitTrace.opto);   % 0/1
nonOptoIdx = optoTrials == 0;
optoIdx    = find(optoTrials == 1);

nNeurons   = numel(allSpk);
SpikesTagged = allSpk;    % initialize output cell array (same size)

for neuronId = 1:nNeurons
    tempSpk = allSpk{neuronId};    % [nTrials x nTime]
    if isempty(tempSpk), continue; end

    if tagged(neuronId)
        % tagged: add opto spikes only in opto trials
        spk_opto_tag = adOpto(tempSpk, optoIdx);      % [nOptoTrials x nTime]

        % start from original spikes, then overwrite opto trials
        spk_combined = tempSpk;
        spk_combined(optoIdx,:) = spk_opto_tag;
    else
        % untagged: keep everything unchanged
        spk_combined = tempSpk;
    end

    % store combined trials (same trial order as original tempSpk)
    SpikesTagged{neuronId} = spk_combined;
end
Spikes.BiPOLES.hit.spks = SpikesTagged;
%% Dataset balencing
allSpk      = Spikes.PSTH.hit.spks;                 % 1 x nNeurons cell
optoTrials  = vertcat(IntanBehaviour.cueHitTrace.opto);  % [nTrials x 1], 0/1

nonOptoIdx  = find(optoTrials == 0);   % numeric indices of non-opto trials
optoIdx     = find(optoTrials == 1);   % numeric indices of opto trials

nNonOpto    = numel(nonOptoIdx);
nOpto       = numel(optoIdx);          % should be 8

% --- choose subsample and upsampling factors ---
targetNonOpto = min(40, nNonOpto);     % e.g. at most 40 non-opto trials
repeatOpto    = max(1,floor(targetNonOpto/length(optoIdx)));   % each opto trial appears x times based on # trials

% random subset of non-opto trial indices
rng(0); % for reproducibility if you like
subNonOptoIdx = randsample(nonOptoIdx, targetNonOpto, false);

% upsampled opto indices by repetition
upsampledOptoIdx = repmat(optoIdx(:), repeatOpto, 1);

nNeurons        = numel(allSpk);
spk_nonOpto_bal = cell(1,nNeurons);  % balanced non-opto (subset)
spk_opto_bal    = cell(1,nNeurons);  % balanced opto (upsampled)

for neuronId = 1:nNeurons
    tempSpk = allSpk{neuronId};      % [nTrials x nTime]
    if isempty(tempSpk), continue; end

    % Subsampled non-opto trials
    spk_nonOpto_bal{neuronId} = tempSpk(subNonOptoIdx, :);

    % Upsampled opto trials (duplicated rows)
    spk_opto_bal{neuronId}    = tempSpk(upsampledOptoIdx, :);
end
nNeurons = numel(spk_nonOpto_bal);
bal_spk  = cell(1, nNeurons);

for neuronId = 1:nNeurons
    spkNon = spk_nonOpto_bal{neuronId};   % [nNonOpto_bal x nTime] or []
    spkOpt = spk_opto_bal{neuronId};      % [nOpto_bal x nTime] or []

    if isempty(spkNon) && isempty(spkOpt)
        bal_spk{neuronId} = [];
    elseif isempty(spkNon)
        bal_spk{neuronId} = spkOpt;
    elseif isempty(spkOpt)
        bal_spk{neuronId} = spkNon;
    else
        % concatenate trials from non-opto and opto along the trial dimension
        bal_spk{neuronId} = [spkNon; spkOpt];   % [ (nNonOpto_bal+nOpto_bal) x nTime ]
    end
end
Spikes.BiPOLES.hit.spks = bal_spk;
%% Analyze neural dynamics of opto and non opto trials
% Note that we concatenate trial conditions as to apply the same models for
% statistical comparison (ie. hit vs miss, hit vs FA, opto vs no opto)
Spikes = makeSpikeGPFA(Spikes);
Spikes.GPFA.HitMiss.dat = [Spikes.GPFA.hit.dat,Spikes.GPFA.miss.dat];
for n = 1:IntanBehaviour.nCueHit%+1:length(Spikes.GPFA.BaselineOpto.dat) %fix trials
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
[neuralDynamics,waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);

%%
optoTrials = vertcat(IntanBehaviour.cueHitTrace.opto);   % [nTrials x 1]

nonOptoIdx = optoTrials == 0;
optoIdx    = optoTrials == 1;

% nonOptoIdx = 1:40;
% optoIdx = 41:76;
% dim x time x trials  ->  time x trials per dim
x = squeeze(neuralDynamics.hitOnly.X(1,:,:));   % [time x trials]
y = squeeze(neuralDynamics.hitOnly.X(2,:,:));
z = squeeze(neuralDynamics.hitOnly.X(3,:,:));

xno = x(:,nonOptoIdx);   yno = y(:,nonOptoIdx);   zno = z(:,nonOptoIdx);
xo  = x(:,optoIdx);      yo  = y(:,optoIdx);      zo  = z(:,optoIdx);

% mean trajectories
mx_no = mean(xno,2); my_no = mean(yno,2); mz_no = mean(zno,2);
mx_o  = mean(xo,2);  my_o  = mean(yo,2);  mz_o  = mean(zo,2);

% colors
col_no = [0.6 0.6 0.6];           % light gray baseline
col_o  = [0 123 167]/255;         % cerulean for opto [web:115]

reactionTime = vertcat(IntanBehaviour.cueHitTrace.reactionTime);
startIdx = 1;
stimIdx  = 75;                    % stimulus bin
rtIdx = floor((1500+mean(reactionTime(nonOptoIdx))*1000)/20);
figure; hold on

% baseline trajectory (thin gray)
plot3(mx_no, my_no, mz_no, 'Color', col_no, 'LineWidth', 2);

% opto trajectory (thicker cerulean)
plot3(mx_o, my_o, mz_o, 'Color', col_o, 'LineWidth', 2.5);

% markers for start and stim non opto trajectory
plot3(mx_no(startIdx), my_no(startIdx), mz_no(startIdx), ...
      'o', 'MarkerSize', 8, 'MarkerFaceColor', col_no, 'MarkerEdgeColor','none');
plot3(mx_no(stimIdx),  my_no(stimIdx),  mz_no(stimIdx),  ...
      'o', 'MarkerSize', 8, 'MarkerFaceColor', col_no, 'MarkerEdgeColor','none');


plot3(mx_no(rtIdx),  my_no(rtIdx),  mz_no(rtIdx),  ...
      'o', 'MarkerSize', 8, 'MarkerFaceColor', col_no, 'MarkerEdgeColor','none');

% markers at start and stim along opto trajectory
plot3(mx_o(startIdx), my_o(startIdx), mz_o(startIdx), ...
      'o', 'MarkerSize', 8, 'MarkerFaceColor', col_o, 'MarkerEdgeColor','none');
plot3(mx_o(stimIdx),  my_o(stimIdx),  mz_o(stimIdx),  ...
      'o', 'MarkerSize', 8, 'MarkerFaceColor', col_o, 'MarkerEdgeColor','none');

rtIdx = floor((1500+mean(reactionTime(optoIdx))*1000)/20);
plot3(mx_o(rtIdx),  my_o(rtIdx),  mz_o(rtIdx),  ...
      'o', 'MarkerSize', 8, 'MarkerFaceColor', col_o, 'MarkerEdgeColor','none');

set(gca,'Box','off','TickDir','out','XColor','k','YColor','k','ZColor','k');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
view(3);
%% Speed
optoTrials = vertcat(IntanBehaviour.cueHitTrace.opto);
nonOptoIdx = optoTrials == 0;
optoIdx    = optoTrials == 1;

% neuralDynamics.hitOnly.speed.speed: dim x time x trials
spd = neuralDynamics.hitOnly.speed.speed;   % assume dim 1 = speed

t = linspace(-1.5,1.5,size(spd,2));   % or construct time vector in s

spd_no = squeeze(spd(1,:,nonOptoIdx));   % [time x nNo]
spd_o  = squeeze(spd(1,:,optoIdx));      % [time x nOp]

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
p_vec = linspace(0.9, 0.2, nPulses);   % low → high probability


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
