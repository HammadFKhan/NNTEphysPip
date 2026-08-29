%% Sparse CC analysis
% Script to generate CCA analysis of neural trajectory during task.

% Following the standard approach in CCA, we identified two sets of loading
% vectors, {wi} and {vi}, termed here as CCA modes, each of which was an
% activity mode within one of the two neural ensembles (that is, with N1
% and N2 elements, respectively). The index i ∈ {1, 2, 3, ..., minimum(N1,
% N2)} denoted the individual modes, which we determined such that the
% projections of the neural activity fluctuations, X and Y, onto wi and vi,
% were maximally correlated between the two trajectories, subject to the
% normalization constraint. Given this normalization condition, the
% quantity ) equals the correlation coefficient of the activity modes,and
% in the two different brain areas. After finding the first CCA mode (i
% =1), we identified successive modes in an iterative manner. Specifically,
% for all previously identified CCA modes we removed the CCA fluctuations
% from X and Y. We applied equation (11) to the residuals and thereby
% identified a set of orthonormal fluctuation modes with correlation
% coefficient values that progressively declined with the index, i. To
% identify the maxima specified by equation (11), we first randomly
% initialized the vectors wi and vi while constraining them to have unity
% length. We then found values of wi and vi that maximized the objective
% function in equation (11) by performing an alternating optimization

% Ebrahimi, S., Lecoq, J., Rumyantsev, O. et al. Emergent reliability in
% sensory cortical coding and inter-area communication. Nature 605, 713–721
% (2022).

% CCA Analysis for eOPN data


addpath(genpath('Main\CCA_utilities'));
% files = dir(fullfile('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\M1M2DualShank\CCA\','*.mat'));
% files = dir(fullfile('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\CCA_cooling','*.mat'));
%files = dir(fullfile('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\M1M2DualShank\CCA\eOPNThal\','*.mat')); %Thal CCA
files = dir(fullfile('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\M1M2DualShank\CCA\eOPNM1\','*.mat')); %Thal CCA

for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    CCAResults(fileNum).filename = files(fileNum).name;
    CCAResults(fileNum).IntanBehaviour = IntanBehaviour;
    % Make M1 and M2 PCA dimensions based on GPFA
    if ~isfield(M1Spikes,'GPFA')
        disp('Running GPFA...')
        M1Spikes = getGPFA(M1Spikes,IntanBehaviour);
        M2Spikes = getGPFA(M2Spikes,IntanBehaviour);
    end
    [M1rh,M1rm,M1rmh,M1rmf] = trajNorm(M1Spikes,IntanBehaviour);
    [M2rh,M2rm,M2rmh,M2rmf] = trajNorm(M2Spikes,IntanBehaviour);

    %     %% Make M1 and M2 based on binned spike data
    %     [M1rh,M1rm,M1rmh,M1rmf] = makeSpikeCCA(M1Spikes);
    %     [M2rh,M2rm,M2rmh,M2rmf] = makeSpikeCCA(M2Spikes);

    % Sparse CCA Analysis
    % Here we take the high dimensional neural trajectory data and perform CCA
    % analysis on it to see what correlations there are from the time varying
    % signals. I chose the top 5 modes based on trajectories that occupy 15
    % latent dimensions.
    % For statistical analysis we build seperate CCA models on subset of trials
    % lets say we only use 80% of the data to check for validity.

    nModes = 5;
    iter = 2; %Number of training rounds
    dataKeep = 0.9; % Percentage we keep for CCA model

    timeLag = NaN;
    timeLagFlag = 0;
    shufFlag = 0;
    IntanBehaviour.parameters.cool = 0;
    IntanBehaviour.parametters.opto = 1;
    %%%  Check for cooling and opto condition
    if IntanBehaviour.parameters.opto == 1

        try
            [IntanBehaviourBaseline,IntanBehaviourOpto, Waves, WavesOpto] = separateOptoTrials(IntanBehaviour,IntanBehaviour.parameters);
        catch
            disp('bad session')
            continue
        end
    baselineId_h = 1:length(IntanBehaviourBaseline.cueHitTrace);
    baselineId_m = 1:length(IntanBehaviourBaseline.cueMissTrace);
    baselineId_FA = 1:length(IntanBehaviourBaseline.MIFATrace);

    eOPNId_h = 1:length(IntanBehaviourOpto.cueHitTrace);
    eOPNId_m = 1:length(IntanBehaviourOpto.cueMissTrace);
    eOPNId_FA = 1:length(IntanBehaviourOpto.MIFATrace);
    if length(baselineId_m) < 1
        baselineId_m = 1:5;
    end
    if length(baselineId_FA) <1
        baselineId_mFA = 1:5;
    end

    eOPNId = length(IntanBehaviourBaseline.cueHitTrace)+1:length(IntanBehaviour.cueHitTrace);

    CCAResults(fileNum).CCABaseline = getCCA(M1rh(:,baselineId_h,:),M1rm(:,baselineId_m,:),M1rmh(:,baselineId_h,:),M1rmf(:,baselineId_FA,:),...
        M2rh(:,baselineId_h,:),M2rm(:,baselineId_m,:),M2rmh(:,baselineId_h,:),M2rmf(:,baselineId_FA,:),iter,nModes,dataKeep,timeLag,shufFlag);

    CCAResults(fileNum).CCAeOPN = getCCA(M1rh(:,eOPNId_h,:),M1rm(:,eOPNId_m,:),M1rmh(:,eOPNId_h,:),M1rmf(:,eOPNId_FA,:),...
        M2rh(:,eOPNId_h,:),M2rm(:,eOPNId_m,:),M2rmh(:,eOPNId_h,:),M2rmf(:,eOPNId_FA,:),iter,nModes,dataKeep,timeLag,shufFlag);
    end

    if IntanBehaviour.parameters.cool
        tempCutoff = -9;
        h = IntanBehaviour.hitTemp>tempCutoff; m = IntanBehaviour.missTemp>tempCutoff;FA = IntanBehaviour.FATemp>tempCutoff;

        CCAResults(fileNum).CCABaseline = getCCA(M1rh(:,h,:),M1rm(:,m,:),M1rmh(:,h,:),M1rmf(:,FA,:),...
            M2rh(:,h,:),M2rm(:,m,:),M2rmh(:,h,:),M2rmf(:,FA,:),iter,nModes,dataKeep,timeLag,shufFlag);

        h = IntanBehaviour.hitTemp<=tempCutoff; m = IntanBehaviour.missTemp<=tempCutoff;FA = IntanBehaviour.FATemp<=tempCutoff;

        CCAResults(fileNum).CCACool = getCCA(M1rh(:,h,:),M1rm(:,m,:),M1rmh(:,h,:),M1rmf(:,FA,:),...
            M2rh(:,h,:),M2rm(:,m,:),M2rmh(:,h,:),M2rmf(:,FA,:),iter,nModes,dataKeep,timeLag,shufFlag);
    end
    if timeLagFlag == 1
        timeLag = -25:5:25;
        CCAResultsLag(fileNum).CCA = getCCA(M1rh,M1rm,M1rmh,M1rmf,M2rh,M2rm,M2rmh,M2rmf,iter,nModes,dataKeep,timeLag,shufFlag);
    end
end
%%
% Control condition where we set the time lag for CCA control. If we set it
% as non negative we let M2 lead M1. If negative then we force M2 to lag
% M1.
timeLag = -25:5:25;
shufFlag = 0;

count = 1;
for n = timeLag
    CCA_timeLags(count).CCA = getCCA(M1rh,M1rm,M1rmh,M1rmf,M2rh,M2rm,M2rmh,M2rmf,iter,nModes,dataKeep,n,shufFlag);
    CCA_timeLags(count).timeLag = n;
    count  = count+1;
end
sessionName = [fpath,'/','CCA_data_thal.mat'];
save(sessionName,"IntanBehaviour","parameters","M1Spikes","M2Spikes","CCA","CCA_timeLags", "fpath","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",

%% Pool CCA responses across sessions

CCAall = [CCAResults];
hitCCATracebaseline = []; missCCATracebaseline = []; MIhitCCATracebaseline = []; MIFACCATracebaseline = [];
hitCCATraceeOPN = []; missCCATraceeOPN = []; MIhitCCATraceeOPN = []; MIFACCATraceeOPN = [];
hitCCATrace = []; missCCATrace = []; MIhitCCATrace = []; MIFACCATrace = [];
for n = 1:length(CCAall)
    try
    dat = arrayfun(@(x) x.rVec(1,:),CCAall(n).CCABaseline.hit,'UniformOutput',false);
    dat = vertcat(dat{:});
    hitCCATracebaseline = vertcat(hitCCATracebaseline,dat);
    hitCCATracebaseline(isnan(hitCCATracebaseline(:,1)),:) = [];
    
%     dat = arrayfun(@(x) x.rVec(1,:),CCAall(n).CCABaseline.miss,'UniformOutput',false);
%     dat = vertcat(dat{:});
%     missCCATracebaseline = vertcat(missCCATracebaseline,dat);
%     missCCATracebaseline(isnan(missCCATracebaseline(:,1)),:) = [];
%     missCCATracebaseline=inpaint_nans(missCCATracebaseline);
% 
% 
%     dat = arrayfun(@(x) x.rVec(1,:),CCAall(n).CCABaseline.MIhit,'UniformOutput',false);
%     dat = vertcat(dat{:});
%     MIhitCCATracebaseline = vertcat(MIhitCCATracebaseline,dat);
%     MIhitCCATracebaseline(isnan(MIhitCCATracebaseline(:,1)),:) = [];
% 
%     dat = arrayfun(@(x) x.rVec(1,:),CCAall(n).CCABaseline.MIFA,'UniformOutput',false);
%     dat = vertcat(dat{:});
%     MIFACCATracebaseline = vertcat(MIFACCATracebaseline,dat);
%     MIFACCATracebaseline(isnan(MIFACCATracebaseline(:,1)),:) = [];

    dat = arrayfun(@(x) x.rVec(1,:),CCAall(n).CCAeOPN.hit,'UniformOutput',false);
    dat = vertcat(dat{:});
    hitCCATraceeOPN = vertcat(hitCCATraceeOPN,dat);
    hitCCATraceeOPN(isnan(hitCCATraceeOPN(:,1)),:) = [];
    
%     dat = arrayfun(@(x) x.rVec(1,:),CCAall(n).CCAeOPN.miss,'UniformOutput',false);
%     dat = vertcat(dat{:});
%     missCCATraceeOPN = vertcat(missCCATraceeOPN,dat);
%     missCCATraceeOPN(isnan(missCCATraceeOPN(:,1)),:) = [];
%     missCCATraceeOPN=inpaint_nans(missCCATraceeOPN);
%     
%     dat = arrayfun(@(x) x.rVec(1,:),CCAall(n).CCAeOPN.MIhit,'UniformOutput',false);
%     dat = vertcat(dat{:});
%     MIhitCCATraceeOPN = vertcat(MIhitCCATraceeOPN,dat);
%     MIhitCCATraceeOPN(isnan(MIhitCCATraceeOPN(:,1)),:) = [];
%     
%     dat = arrayfun(@(x) x.rVec(1,:),CCAall(n).CCAeOPN.MIFA,'UniformOutput',false);
%     dat = vertcat(dat{:});
%     MIFACCATraceeOPNe = vertcat(MIFACCATraceeOPN,dat);
%     MIFACCATraceeOPNe(isnan(MIFACCATraceeOPNe(:,1)),:) = [];

    catch 
        disp(['CCA not computed for session ' num2str(n)])
        continue
    end
end

%% ---- Preprocess CCA traces (baseline + smoothing only) ----
baselineIdx = 1:70;
smoothWin   = 3;

missCCATrace=inpaint_nans(missCCATrace);
hitCCA_proc_baseline   = preprocessCCA(hitCCATracebaseline,    baselineIdx, smoothWin);
% missCCA_proc_baseline  = preprocessCCA(missCCATracebaseline,   baselineIdx, smoothWin);
% MIHitCCA_proc_baseline = preprocessCCA(MIhitCCATracebaseline,  baselineIdx, smoothWin);
% FACCA_proc_baseline    = preprocessCCA(MIFACCATracebaseline,   baselineIdx, smoothWin);

hitCCA_proc_eOPN   = preprocessCCA(hitCCATraceeOPN,    baselineIdx, smoothWin);
% missCCA_proc_eOPN  = preprocessCCA(missCCATraceeOPN,   baselineIdx, smoothWin);
% MIHitCCA_proc_eOPN = preprocessCCA(MIhitCCATraceeOPN,  baselineIdx, smoothWin);
% FACCA_proc_eOPN    = preprocessCCA(MIFACCATraceeOPN,   baselineIdx, smoothWin);

t = linspace(-1.5,1.5,size(hitCCA_proc_baseline,2));
stimIdx = 75;

plotCCAeOPNconditions(hitCCA_proc_eOPN,hitCCA_proc_baseline, t, stimIdx);
title('Baseline vs eOPN')
xlim([-0.5 1.5])
axis square
%plotCCAconditions(hitCCA_proc, missCCA_proc, MIHitCCA_proc, FACCA_proc, t, stimIdx);


%%
% hitCCA_proc, missCCA_proc, MIHitCCA_proc, FACCA_proc: [nTrials x nTime]
% Define pre- and post-stim windows (indices into time axis)
preIdx  = 1:75;        % pre-stim (adjust as needed)
postIdx = 75:150;   % post-stim

% Pre/post difference per trial: |post mean - pre mean|
hitCCApost = mean(hitCCA_proc_baseline(:,postIdx),2);
[~,id] = maxk(hitCCApost,5);
hitCCA   = abs(mean(hitCCA_proc_baseline(:,postIdx), 2, 'omitnan') - ...
               mean(hitCCA_proc_baseline(:,preIdx),  2, 'omitnan'));

eOPNCCA = abs(mean(hitCCA_proc_eOPN(:,postIdx), 2, 'omitnan') - ...
               mean(hitCCA_proc_eOPN(:,preIdx),  2, 'omitnan'));


% hitCCA, missCCA, MIhitCCA, MIFACCA are column vectors

[~, p_hit_miss] = ttest(hitCCA, eOPNCCA);   % paired t-test
disp(p_hit_miss)

figure;
plotPairedBarScatter(hitCCA, eOPNCCA, ...
                     'Hit', 'eOPN', 'CC Coefficient', p_hit_miss);



%% LOCAL FUNCTIONS
function plotPairedBarScatter(dataA, dataB, labelA, labelB, yLabelStr,pVal)
% dataA, dataB: column vectors or same-length row vectors (paired)
% labelA, labelB: strings for x‑axis labels
% yLabelStr: y‑axis label
% pVal (optional): p‑value to display above the bars

dataA = dataA(:);
dataB = dataB(:);
assert(numel(dataA) == numel(dataB), 'dataA and dataB must have same length');

n = numel(dataA);
xA = ones(n,1);
xB = 2*ones(n,1);

meanA = mean(dataA,'omitnan');
meanB = mean(dataB,'omitnan');

hold on;

% bar colors (light gray bars)
barWidth = 0.5;
bar(1, meanA, barWidth, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none');
bar(2, meanB, barWidth, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none');

% scatter + connecting lines
for i = 1:n
    plot([xA(i) xB(i)], [dataA(i) dataB(i)], '-', 'Color', [0.6 0.6 0.6]); % line
end
scatter(xA, dataA, 35, [0.10 0.45 0.85], 'filled');  % blue-ish
scatter(xB, dataB, 35, [0.4 0.4 0.4], 'filled');     % gray

% axis formatting
xlim([0.5 2.5]);
xticks([1 2]);
xticklabels({labelA, labelB});
ylabel(yLabelStr);
box off;
set(gca,'TickDir','out','FontSize',12);

% y‑limits with a bit of headroom
yl = ylim;
ylim([0 max(yl(2), max([dataA; dataB])*1.1)]);

% optional p‑value text
if nargin >= 6 && ~isempty(pVal)
    yText = ylim;
    yText = yText(2) * 1.02;
    text(1.5, yText, sprintf('%.1g', pVal), ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'FontSize',10);
end
end

function CCA_proc = preprocessCCA(CCA_in, baselineIdx, smoothWin)
% CCA_in:   [nEpochs x nTime] CCA traces
% baselineIdx: indices used for baseline (e.g., 1:50)
% smoothWin:   moving-average window (samples), e.g., 3

CCA_in = double(CCA_in);             % ensure numeric
[nEpochs, nTime] = size(CCA_in);

%% 1) Baseline correction per epoch
CCA_bc = CCA_in;
for e = 1:nEpochs
    base = mean(CCA_in(e, baselineIdx), 'omitnan');
    CCA_bc(e,:) = CCA_in(e,:) - base;
end

%% 2) Smoothing (zero-phase moving average)
if smoothWin > 1
    k = ones(1, smoothWin) / smoothWin;
    CCA_s = zeros(size(CCA_bc));
    for e = 1:nEpochs
        CCA_s(e,:) = filtfilt(k, 1, CCA_bc(e,:));
    end
else
    CCA_s = CCA_bc;
end

%% 3) Normalize each epoch to [0,1]
CCA_proc = zeros(size(CCA_s));
for e = 1:nEpochs
    x = CCA_s(e,:);
    xMin = min(x, [], 'omitnan');
    xMax = max(x, [], 'omitnan');
    if xMax > xMin
        CCA_proc(e,:) = (x - xMin) / (xMax - xMin);
    else
        CCA_proc(e,:) = zeros(1, nTime);   % flat trace
    end
end
end

function plotCCAconditions(hitCCA_proc, missCCA_proc, MIHitCCA_proc, FACCA_proc, t, stimIdx)

% Colors
colHit   = [0.2 0.4 0.9];
colMiss  = [0.4 0.4 0.4];
colMIHit = [0.2 0.4 0.9];
colFA    = [0.8 0.2 0.2];

%% ---- Compute mean and SEM ----
meanHit    = mean(hitCCA_proc,1,'omitnan');
seHit      = std(hitCCA_proc,0,1,'omitnan') ./ sqrt(size(hitCCA_proc,1));

meanMiss   = mean(missCCA_proc,1,'omitnan');
seMiss     = std(missCCA_proc,0,1,'omitnan') ./ sqrt(size(missCCA_proc,1));

meanMIHit  = mean(MIHitCCA_proc,1,'omitnan');
seMIHit    = std(MIHitCCA_proc,0,1,'omitnan') ./ sqrt(size(MIHitCCA_proc,1));

meanFA     = mean(FACCA_proc,1,'omitnan');
seFA       = std(FACCA_proc,0,1,'omitnan') ./ sqrt(size(FACCA_proc,1));

%% ---- Hit vs Miss (baseline-aligned) ----
preIdx = 1:75;   % pre‑stim window for condition-level baseline

baseHit  = mean(meanHit(preIdx),  'omitnan');
baseMiss = mean(meanMiss(preIdx), 'omitnan');
commonBase = (baseHit + baseMiss) / 2;

meanHit_shift  = meanHit  - baseHit  + commonBase;
meanMiss_shift = meanMiss - baseMiss + commonBase;
seHit_shift    = seHit;
seMiss_shift   = seMiss;

figure; hold on;
% Hit
fill([t fliplr(t)], [meanHit_shift+seHit_shift ...
                     fliplr(meanHit_shift-seHit_shift)], ...
     colHit, 'FaceAlpha',0.25, 'EdgeColor','none');
plot(t, meanHit_shift, 'Color', colHit, 'LineWidth', 2);
% Miss
fill([t fliplr(t)], [meanMiss_shift+seMiss_shift ...
                     fliplr(meanMiss_shift-seMiss_shift)], ...
     colMiss, 'FaceAlpha',0.25, 'EdgeColor','none');
plot(t, meanMiss_shift, 'Color', colMiss, 'LineWidth', 2);

xline(stimIdx, '--', 'Color', [0.6 0 0], 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Normalized CCA (0–1, aligned baseline)');
legend({'Hit \pm SEM','Hit mean','Miss \pm SEM','Miss mean'}, 'Location','best');
title('CCA traces: Hit vs Miss (baseline-aligned)');
box off; set(gca,'TickDir','out','FontSize',12);
xlim([-1.0 max(t)]);

%% ---- MI Hit vs MI False Alarm (baseline-aligned) ----
baseMIHit = mean(meanMIHit(preIdx), 'omitnan');
baseFA    = mean(meanFA(preIdx),    'omitnan');
commonBase_MI = (baseMIHit + baseFA) / 2;

meanMIHit_shift = meanMIHit - baseMIHit + commonBase_MI;
meanFA_shift    = meanFA    - baseFA    + commonBase_MI;
seMIHit_shift   = seMIHit;
seFA_shift      = seFA;

figure; hold on;
% MI Hit
fill([t fliplr(t)], [meanMIHit_shift+seMIHit_shift ...
                     fliplr(meanMIHit_shift-seMIHit_shift)], ...
     colMIHit, 'FaceAlpha',0.25, 'EdgeColor','none');
plot(t, meanMIHit_shift, 'Color', colMIHit, 'LineWidth', 2);
% MI FA
fill([t fliplr(t)], [meanFA_shift+seFA_shift ...
                     fliplr(meanFA_shift-seFA_shift)], ...
     colFA, 'FaceAlpha',0.25, 'EdgeColor','none');
plot(t, meanFA_shift, 'Color', colFA, 'LineWidth', 2);

xline(stimIdx, '--', 'Color', [0.6 0 0], 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Normalized CCA (0–1, aligned baseline)');
legend({'MI Hit \pm SEM','MI Hit mean','MI FA \pm SEM','MI FA mean'}, ...
       'Location','best');
title('CCA traces: MI Hit vs MI False Alarm (baseline-aligned)');
box off; set(gca,'TickDir','out','FontSize',12);
xlim([-1.0 max(t)]);
end

function plotCCAeOPNconditions(hitCCA_proc, eOPNCCA_proc, t, stimIdx)

% Colors
colHit   = [0.2 0.4 0.9];
coleOPN  = [217/255 83/255 25/255];


%% ---- Compute mean and SEM ----
meanHit    = mean(hitCCA_proc,1,'omitnan');
seHit      = std(hitCCA_proc,0,1,'omitnan') ./ sqrt(size(hitCCA_proc,1)*2);

meaneOPN   = mean(eOPNCCA_proc,1,'omitnan');
seeOPN     = std(eOPNCCA_proc,0,1,'omitnan') ./ sqrt(size(eOPNCCA_proc,1)*2);



%% ---- Hit vs Miss (baseline-aligned) ----
preIdx = 1:70;   % pre‑stim window for condition-level baseline

baseHit  = mean(meanHit(preIdx),  'omitnan');
baseMiss = mean(meaneOPN(preIdx), 'omitnan');
commonBase = (baseHit + baseMiss) / 2;

meanHit_shift  = meanHit  - baseHit  + commonBase;
meaneOPN_shift = meaneOPN - baseMiss + commonBase;
seHit_shift    = seHit;
seMiss_shift   = seeOPN;

figure; hold on;
% Hit
fill([t fliplr(t)], [meanHit_shift+seHit_shift ...
                     fliplr(meanHit_shift-seHit_shift)], ...
     colHit, 'FaceAlpha',0.25, 'EdgeColor','none');
plot(t, meanHit_shift, 'Color', colHit, 'LineWidth', 2);
% Miss
fill([t fliplr(t)], [meaneOPN_shift+seMiss_shift ...
                     fliplr(meaneOPN_shift-seMiss_shift)], ...
     coleOPN, 'FaceAlpha',0.25, 'EdgeColor','none');
plot(t, meaneOPN_shift, 'Color', coleOPN, 'LineWidth', 2);

xline(stimIdx, '--', 'Color', [0.6 0 0], 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Normalized CCA (0–1, aligned baseline)');
legend({'Hit \pm SEM','Hit mean','eOPN \pm SEM','eOPN mean'}, 'Location','best');
title('CCA traces: Hit vs Miss (baseline-aligned)');
box off; set(gca,'TickDir','out','FontSize',12);
xlim([-1.0 max(t)]);
end

%% LOCAL FUNCTIONS
function r = meanTraj(X,trials,components)
r = X(:,:,trials);
r = permute(r,[1 3 2]);
end

function [rh,rm,rmh,rmf] = trajNorm(Spikes,Behaviour)
X = arrayfun(@(x) vertcat(x.xorth),Spikes.GPFA.seqTrainHitMiss,'UniformOutput',false);
X = horzcat(X{:});
neuralTrajHitMiss = reshape(X,size(X,1),Spikes.GPFA.seqTrainHitMiss(1).T,[]);

X = neuralTrajHitMiss;
hittrials = 1:length(Behaviour.cueHitTrace);
misstrials = length(Behaviour.cueHitTrace)+1:size(X,3);
rh = meanTraj(X,hittrials,6); %trajectory variable and predefined conditional trial indexes
rm = meanTraj(X,misstrials,6); %trajectory variable and predefined conditional trial indexes

X = arrayfun(@(x) vertcat(x.xorth),Spikes.GPFA.seqTrainMIHitFA,'UniformOutput',false);
X = horzcat(X{:});
neuralTrajMIHitFA = reshape(X,size(X,1),Spikes.GPFA.seqTrainMIHitFA(1).T,[]);

X = neuralTrajMIHitFA;

hittrials = 1:length(Behaviour.MIHitTrace);
FAtrials = length(Behaviour.MIHitTrace)+1:size(X,3);

rmh = meanTraj(X,hittrials,6);
rmf = meanTraj(X,FAtrials,6);

end

function CCA = getCCA(M1rh,M1rm,M1rmh,M1rmf,M2rh,M2rm,M2rmh,M2rmf,iter,nModes,dataKeep,timeLagArray,shuf)
% Initialize Structure
CCA = struct();
% CCA.hit.wxMat = [];     CCA.miss.wxMat = [];       CCA.MIhit.wxMat = [];    CCA.MIFA.wxMat = [];
% CCA.hit.wyMat = [];     CCA.miss.wyMat = [];       CCA.MIhit.wyMat = [];    CCA.MIFA.wyMat = [];
% CCA.hit.rVec = [];      CCA.miss.rVec = [];        CCA.MIhit.rVec = [];    CCA.MIFA.rVec = [];

for nn = 1:iter
    datIdx = randperm(size(M1rh,2));
    datIdx = datIdx(1:floor(length(datIdx)*dataKeep));
    xh = squeeze(M1rh(:,sort(datIdx),:));
    yh = squeeze(M2rh(:,sort(datIdx),:));

    datIdx = randperm(size(M1rm,2));
    datIdx = datIdx(1:floor(length(datIdx)*dataKeep));
    xm = squeeze(M1rm(:,sort(datIdx),:));
    ym = squeeze(M2rm(:,sort(datIdx),:));

    datIdx = randperm(size(M1rmh,2));
    datIdx = datIdx(1:floor(length(datIdx)*dataKeep));
    xmh = squeeze(M1rmh(:,sort(datIdx),:));
    ymh = squeeze(M2rmh(:,sort(datIdx),:));

    datIdx = randperm(size(M1rmf,2));
    datIdx = datIdx(1:floor(length(datIdx)*dataKeep));
    xmf = squeeze(M1rmf(:,sort(datIdx),:));
    ymf = squeeze(M2rmf(:,sort(datIdx),:));


    for n = 1:size(M1rh,3)
        if ~isnan(timeLagArray)
            for timeLagId = 1:length(timeLagArray)
                timeLag = timeLagArray(timeLagId);
                if timeLag>0 %Checks if we want to do time lags
                    %                     disp('Case 1: M2 leads')
                    Xh = squeeze(xh(:,:,n));
                    Xm = squeeze(xm(:,:,n));
                    Xmh = squeeze(xmh(:,:,n));
                    Xmf = squeeze(xmf(:,:,n));
                    if (n+timeLag)<size(M1rh,3)
                        Yh = squeeze(yh(:,:,n+timeLag)); %lead M2 by a certain amount
                        Ym = squeeze(ym(:,:,n+timeLag));
                        Ymh = squeeze(ymh(:,:,n+timeLag));
                        Ymf = squeeze(ymf(:,:,n+timeLag));
                    else
                        Yh = squeeze(yh(:,:,abs(size(M1rh,3)-n+timeLag))); %unless we reach the end of the timepoints
                        Ym = squeeze(ym(:,:,abs(size(M1rh,3)-n+timeLag)));
                        Ymh = squeeze(ymh(:,:,abs(size(M1rh,3)-n+timeLag)));
                        Ymf = squeeze(ymf(:,:,abs(size(M1rh,3)-n+timeLag)));
                    end

                elseif timeLag<0
                    %                     disp('Case 2: M1 leads')
                    if (n+abs(timeLag))<size(M1rh,3)
                        Xh = squeeze(xh(:,:,n+abs(timeLag))); %lead M1 by a certain amount
                        Xm = squeeze(xm(:,:,n+abs(timeLag)));
                        Xmh = squeeze(xmh(:,:,n+abs(timeLag)));
                        Xmf = squeeze(xmf(:,:,n+abs(timeLag)));
                    else
                        Xh = squeeze(xh(:,:,abs(size(M1rh,3)-n+abs(timeLag)))); %unless we reach the end of the timepoints
                        Xm = squeeze(xm(:,:,abs(size(M1rh,3)-n+abs(timeLag))));
                        Xmh = squeeze(xmh(:,:,abs(size(M1rh,3)-n+abs(timeLag))));
                        Xmf = squeeze(xmf(:,:,abs(size(M1rh,3)-n+abs(timeLag))));
                    end
                    Yh = squeeze(yh(:,:,n));
                    Ym = squeeze(ym(:,:,n));
                    Ymh = squeeze(ymh(:,:,n));
                    Ymf = squeeze(ymf(:,:,n));
                end
                %                 disp(['Time Lag ' num2str(timeLag) 'ms'])
                try
                    [CCA.hit(nn).timelag(timeLagId).wxMat(:,:,n),CCA.hit(nn).timelag(timeLagId).wyMat(:,:,n),CCA.hit(nn).timelag(timeLagId).rVec(:,n)]=SparseCCA(Xh',Yh',2,2,1,nModes);
                    [CCA.miss(nn).timelag(timeLagId).wxMat(:,:,n),CCA.miss(nn).timelag(timeLagId).wyMat(:,:,n),CCA.miss(nn).timelag(timeLagId).rVec(:,n)]=SparseCCA(Xm',Ym',2,2,1,nModes);
                    [CCA.MIhit(nn).timelag(timeLagId).wxMat(:,:,n),CCA.MIhit(nn).timelag(timeLagId).wyMat(:,:,n),CCA.MIhit(nn).timelag(timeLagId).rVec(:,n)]=SparseCCA(Xmh',Ymh',2,2,1,nModes);
                    [CCA.MIFA(nn).timelag(timeLagId).wxMat(:,:,n),CCA.MIFA(nn).timelag(timeLagId).wyMat(:,:,n),CCA.MIFA(nn).timelag(timeLagId).rVec(:,n)]=SparseCCA(Xmf',Ymf',2,2,1,nModes);
                    CCA.hit(nn).timelag(timeLagId).timeLag = timeLag;
                    CCA.miss(nn).timelag(timeLagId).timeLag = timeLag;
                    CCA.MIhit(nn).timelag(timeLagId).timeLag = timeLag;
                    CCA.MIFA(nn).timelag(timeLagId).timeLag = timeLag;
                catch
                    disp('error calculating CCA...')
                    continue
                end
                disp(['Timestep ' num2str(n) ' on iteration ' num2str(nn) '...'])
            end
        else
            % Ensure each array is 3D (dims x trials x time)
            xh  = ensure3D(xh);
            xm  = ensure3D(xm);
            xmh = ensure3D(xmh);
            xmf = ensure3D(xmf);
            yh  = ensure3D(yh);
            ym  = ensure3D(ym);
            ymh = ensure3D(ymh);
            ymf = ensure3D(ymf);

            % Now safely index with n
            Xh  = squeeze(xh(:,:,n));
            Xm  = squeeze(xm(:,:,n));
            Xmh = squeeze(xmh(:,:,n));
            Xmf = squeeze(xmf(:,:,n));
            Yh  = squeeze(yh(:,:,n));
            Ym  = squeeze(ym(:,:,n));
            Ymh = squeeze(ymh(:,:,n));
            Ymf = squeeze(ymf(:,:,n));
            try
                [CCA.hit(nn).wxMat(:,:,n),CCA.hit(nn).wyMat(:,:,n),CCA.hit(nn).rVec(:,n)]=SparseCCA(Xh',Yh',2,2,1,nModes);
                [CCA.miss(nn).wxMat(:,:,n),CCA.miss(nn).wyMat(:,:,n),CCA.miss(nn).rVec(:,n)]=SparseCCA(Xm',Ym',2,2,1,nModes);
                [CCA.MIhit(nn).wxMat(:,:,n),CCA.MIhit(nn).wyMat(:,:,n),CCA.MIhit(nn).rVec(:,n)]=SparseCCA(Xmh',Ymh',2,2,1,nModes);
                [CCA.MIFA(nn).wxMat(:,:,n),CCA.MIFA(nn).wyMat(:,:,n),CCA.MIFA(nn).rVec(:,n)]=SparseCCA(Xmf',Ymf',2,2,1,nModes);
            catch
                disp('error calculating CCA...')
                continue
            end
            disp(['Timestep ' num2str(n) ' on iteration ' num2str(nn) '...'])
        end
        %         if shuf
        %             Xh = squeeze(xh(:,:,randperm(size(M1rh,3),1)));
        %             Xm = squeeze(xm(:,:,randperm(size(M1rm,3),1)));
        %             Xmh = squeeze(xmh(:,:,randperm(size(M1rmh,3),1)));
        %             Xmf = squeeze(xmf(:,:,randperm(size(M1rmf,3),1)));
        %
        %             Yh = squeeze(yh(:,:,randperm(size(M1rm,3),1)));
        %             Ym = squeeze(ym(:,:,randperm(size(M1rm,3),1)));
        %             Ymh = squeeze(ymh(:,:,randperm(size(M1rm,3),1)));
        %             Ymf = squeeze(ymf(:,:,randperm(size(M1rm,3),1)));
        %         end

    end
end
% Save parameters
CCA.params.iter = iter;
CCA.params.nModes = nModes;
CCA.params.dataKeep = dataKeep;
CCA.params.timeLag = timeLagArray;
CCA.params.shufFlag = shuf;
end

function [rh,rm,rmh,rmf] = makeSpikeCCA(Spikes)
rh = horzcat(Spikes.rawPSTH.hit.spks{:});
rh = reshape(rh,size(Spikes.rawPSTH.hit.spks{1},1),size(Spikes.rawPSTH.hit.spks{1},2),[]);
rh = permute(rh,[3 1 2]);
kernal = [0 0 0;];
temp = [];
fprintf('Cleaning up hit spikes...\n')
for n = 1:size(rh,1)
    dat  = squeeze(rh(n,:,:));
    dat(:,1500:end) = conv2(dat(:,1500:end),kernal,'same');
    dat = smoothdata(dat,2,'gaussian',15);
    temp(n,:,:) = dat;
end
rh = temp;

fprintf('Cleaning up miss spikes...\n')
rm = horzcat(Spikes.rawPSTH.miss.spks{:});
rm = reshape(rm,size(Spikes.rawPSTH.miss.spks{1},1),size(Spikes.rawPSTH.miss.spks{1},2),[]);
rm = permute(rm,[3 1 2]);
temp = [];
for n = 1:size(rm,1)
    dat  = squeeze(rm(n,:,:));
    dat(:,1500:end) = conv2(dat(:,1500:end),kernal,'same');
    dat = smoothdata(dat,2,'gaussian',25);
    temp(n,:,:) = dat;
end
rm = temp;

fprintf('Cleaning up MI spikes...\n')
rmh = horzcat(Spikes.rawPSTH.MIHit.spks{:});
rmh = reshape(rmh,size(Spikes.rawPSTH.MIHit.spks{1},1),size(Spikes.rawPSTH.MIHit.spks{1},2),[]);
rmh = permute(rmh,[3 1 2]);
temp = [];
for n = 1:size(rmh,1)
    dat  = squeeze(rmh(n,:,:));
    dat(:,1500:end) = conv2(dat(:,1500:end),kernal,'same');
    dat = smoothdata(dat,2,'gaussian',25);
    temp(n,:,:) = dat;
end
rmh = temp;

fprintf('Cleaning up FA spikes...\n')
rmf = horzcat(Spikes.rawPSTH.MIFA.spks{:});
rmf = reshape(rmf,size(Spikes.rawPSTH.MIFA.spks{1},1),size(Spikes.rawPSTH.MIFA.spks{1},2),[]);
rmf = permute(rmf,[3 1 2]);
temp = [];
for n = 1:size(rmf,1)
    dat  = squeeze(rmf(n,:,:));
    dat(:,1500:end) = conv2(dat(:,1500:end),kernal,'same');
    dat = smoothdata(dat,2,'gaussian',25);
    temp(n,:,:) = dat;
end
rmf = temp;
end

function A = ensure3D(A)
sz = size(A);
if ndims(A) == 2
    % Interpret as (dims x time) → make (dims x 1 x time)
    A = reshape(A, sz(1), 1, sz(2));
end
end

function Spikes = getGPFA(Spikes,IntanBehaviour)
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
try
    [Spikes.GPFA.resultHit,Spikes.GPFA.seqTrainHit] = gpfaAnalysis(Spikes.GPFA.hit.dat,1); %Run index
    [Spikes.GPFA.resultMiss,Spikes.GPFA.seqTrainMiss] = gpfaAnalysis(Spikes.GPFA.miss.dat,2); %Run index
    [Spikes.GPFA.resultMIHit,Spikes.GPFA.seqTrainMIHit] = gpfaAnalysis(Spikes.GPFA.MIHit.dat,3); %Run index
    [Spikes.GPFA.resultMIFA,Spikes.GPFA.seqTrainMIFA] = gpfaAnalysis(Spikes.GPFA.MIFA.dat,4); %Run index
    [Spikes.GPFA.resultHitMiss,Spikes.GPFA.seqTrainHitMiss] = gpfaAnalysis(Spikes.GPFA.HitMiss.dat,5); %Run index
    [Spikes.GPFA.resultMIHitFA,Spikes.GPFA.seqTrainMIHitFA] = gpfaAnalysis(Spikes.GPFA.MIHitFA.dat,6); %Run index
catch ME
    disp('Error running GPFA, skipping....')
end
close all
end