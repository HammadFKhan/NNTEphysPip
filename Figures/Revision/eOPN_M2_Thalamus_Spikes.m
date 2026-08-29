%% Generate figure related to pooled eOPN inactivation experiments in primary motor cortex and primary motor thalamus 
% Combining eOPN data together
M2eOPN = struct();
dirFiles = "Y:\Hammad\Ephys\LeverTask\LeverTaskRebuttal\eOPN3\eOPN_dual_thalamus\eOPN_Th_files.xlsx";

% Read the data as a table (recommended)
T = readtable(dirFiles);

colNames = T.Properties.VariableNames;  % column names as cell array of strings
Spks_to_analyze = T(T.Spks == 1 & T.Region == "M1", :); % Change to M1 or Th depending on region
 %%
totalSpikes.nonOptoHit = [];
totalSpikes.optoHit = [];
for fileNum = 1:size(Spks_to_analyze,1)
    disp(['File number: ' num2str(fileNum)])
    fileName = [char(Spks_to_analyze.FilePath(fileNum)) '\' char(Spks_to_analyze.IntanFileName(fileNum)) '\UCLA_chanmap_64F2\Spikes.mat'];
    load(fileName)
%     if ~isfield(Spikes,'GPFA')
%         try
%         Spikes = makeSpikeGPFA(Spikes);
%         Spikes.GPFA.HitMiss.dat = [Spikes.GPFA.hit.dat,Spikes.GPFA.miss.dat];
%         for n = 1:IntanBehaviour.nCueHit%+1:length(Spikes.GPFA.BaselineOpto.dat) %fix trials
%             Spikes.GPFA.HitMiss.dat(n).trialId = n;
%         end
%         Spikes.GPFA.MIHitFA.dat = [Spikes.GPFA.MIHit.dat,Spikes.GPFA.MIFA.dat];
%         for n = length(IntanBehaviour.MIHitTrace)+1:length(Spikes.GPFA.MIHitFA.dat) %fix trials
%             Spikes.GPFA.MIHitFA.dat(n).trialId = n;
%         end
%         %%%
%         addpath(genpath('C:\Users\khan332\Documents\GitHub\NeuralTraj'));
%         addpath(genpath('mat_results'));
%         if exist('mat_results','dir'),rmdir('mat_results','s'),end
%         [Spikes.GPFA.resultHit,Spikes.GPFA.seqTrainHit] = gpfaAnalysis(Spikes.GPFA.hit.dat,1); %Run index
%         [Spikes.GPFA.resultMiss,Spikes.GPFA.seqTrainMiss] = gpfaAnalysis(Spikes.GPFA.miss.dat,2); %Run index
%         [Spikes.GPFA.resultMIHit,Spikes.GPFA.seqTrainMIHit] = gpfaAnalysis(Spikes.GPFA.MIHit.dat,3); %Run index
%         [Spikes.GPFA.resultMIFA,Spikes.GPFA.seqTrainMIFA] = gpfaAnalysis(Spikes.GPFA.MIFA.dat,4); %Run index
%         [Spikes.GPFA.resultHitMiss,Spikes.GPFA.seqTrainHitMiss] = gpfaAnalysis(Spikes.GPFA.HitMiss.dat,5); %Run index
%         [Spikes.GPFA.resultMIHitFA,Spikes.GPFA.seqTrainMIHitFA] = gpfaAnalysis(Spikes.GPFA.MIHitFA.dat,6); %Run index
%         close all
%         catch ME
%             disp('Error calculating GPFA...')
%             continue
%         end
%     end
%     
%     [M1eOPN(fileNum).neuralDynamics,waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);
    M2eOPN(fileNum).IntanBehaviour = IntanBehaviour;

    try
        [IntanBehaviourBaseline,IntanBehaviourOpto, Waves, WavesOpto] = separateOptoTrials(IntanBehaviour,IntanBehaviour.parameters);
        baselineId = 1:length(IntanBehaviourBaseline.cueHitTrace);
        eOPNId = length(IntanBehaviourBaseline.cueHitTrace)+1:length(IntanBehaviour.cueHitTrace);
        assert(length(eOPNId)==length(IntanBehaviourOpto.cueHitTrace))
        % Since the sqEntropy is accessing each spike structure we need to seperate the trials
        baselineSpikes = Spikes;
        eOPNSpikes = Spikes;
        baselineSpikes.PSTH.hit.spks = cellfun(@(x) x(baselineId,:),baselineSpikes.PSTH.hit.spks,'UniformOutput',false);
        eOPNSpikes.PSTH.hit.spks = cellfun(@(x) x(eOPNId,:),eOPNSpikes.PSTH.hit.spks,'UniformOutput',false);
        Spikes.baselineSpikes = baselineSpikes;
        Spikes.eOPNSpikes = eOPNSpikes;

        M2eOPN(fileNum).Spikes = Spikes;
    catch
        disp('error on eopn ')
        continue
    end

    trials = baselineSpikes.PSTH.hit.spks;
    
    output = make_nice_mean_raster(trials,20,0);
    totalSpikes.nonOptoHit = [totalSpikes.nonOptoHit;output];
    

    trials = eOPNSpikes.PSTH.hit.spks;
    output = make_nice_mean_raster(trials,20,0);
    totalSpikes.optoHit = [totalSpikes.optoHit;output];

end
%% Total spikes
spikeRate = smoothdata(totalSpikes.nonOptoHit,2,'gaussian',50);
[hitnormSpk,hittimIdx,hitspkIdx] = spknorm(spikeRate);
f = figure,subplot(131)
plotSpkSeq(hitnormSpk(:,1000:end))
title('Hit')
colormap(flip(gray))
set(gca,'fontsize',16)

spikeRate = smoothdata(totalSpikes.optoHit,2,'gaussian',50);
[hitnormSpk,hittimIdx,~] = spknorm(spikeRate);
f = figure,subplot(131)

plotSpkSeq(hitnormSpk(:,1000:end),hitspkIdx)
title('Hit')
blues = slanCM('Oranges')
colormap(blues)
set(gca,'fontsize',16)
%% Analyze baseline and eopn spikes
% Analyze baseline vs eopn spikes
dynamics = M2eOPN;

win1    = 1500;    % ms indices in original 1‑ms bins
win2    = 1700;
binSize = 50;      % ms per analysis bin

FR_baseline_all = [];   % aggregated across sessions
FR_eopn_all     = [];

for n = 1:numel(dynamics)

    baselineSpks = dynamics(n).Spikes.baselineSpikes.PSTH.hit.spks;   % 1 x nUnits cell
    eopnSpks     = dynamics(n).Spikes.eOPNSpikes.PSTH.hit.spks;      % 1 x nUnits cell

    nUnits = numel(baselineSpks);
    fr_baseline = nan(1,nUnits);
    fr_eopn     = nan(1,nUnits);

    for u = 1:nUnits
        % each cell: [nTrials x nTime]
        B = baselineSpks{u}(:,win1:win2);
        E = eopnSpks{u}(:,win1:win2);

        nTime = size(B,2);
        edges = 1:binSize:nTime+1;
        nBins = numel(edges)-1;

        cntB = zeros(size(B,1),nBins);
        cntE = zeros(size(E,1),nBins);

        for b = 1:nBins
            idx = edges(b):edges(b+1)-1;
            cntB(:,b) = sum(B(:,idx),2);
            cntE(:,b) = sum(E(:,idx),2);
        end

        % spikes per bin -> Hz [web:19]
        rateB = cntB * (1000/binSize);
        rateE = cntE * (1000/binSize);
        fr_baseline(u) = mean(rateB(:));
        fr_eopn(u)     = mean(rateE(:));
    end

    % append this session’s units to the global arrays [web:54][web:50]
    FR_baseline_all = [FR_baseline_all, fr_baseline];
    FR_eopn_all     = [FR_eopn_all,     fr_eopn];
end
% plot it out
f = figure; hold on
f.Position = [1000         250         560         420];
scatter(FR_baseline_all, FR_eopn_all, 25, 'k', 'filled');
lims = [0 max([FR_baseline_all FR_eopn_all])*1.05];
plot(lims, lims, '--', 'Color', [0.3 0.8 0.6], 'LineWidth', 1.5);
axis square; xlim(lims); ylim(lims);
xlabel('spikes s^{-1} (control)');
ylabel('spikes s^{-1} (light)');


axis square                          % equal x/y scale [web:11]
xlim(lims); ylim(lims);              % same limits on both axes [web:11][web:17]
xlabel('spikes s^{-1} (control)','FontSize',10);  % [web:12]
ylabel('spikes s^{-1} (light)','FontSize',10,...
       'Color',[0.6 0 0.6]);                                             % magenta-ish y label [web:12]

set(gca,'Box','off','TickDir','out','FontSize',9);
title('M2 to M1 eOPN3');

% Wilcoxon signed-rank test (paired, two-sided)
[p,h,stats] = signrank(FR_baseline_all, FR_eopn_all);   % [web:44]

fprintf('Wilcoxon signed-rank test: p = %.3g, z = %.3f\n', p, stats.zval);
txt = sprintf('Wilcoxon signed-rank: p = %.3g', p);
 
txt2 = sprintf('Baseline FR: %.3g, eOPN FR: %.3g\n ', mean(fr_baseline),mean(fr_eopn));
text(0.05*max(xlim), 0.9*max(ylim), txt, 'FontSize', 9);
text(0.05*max(xlim), 0.8*max(ylim), txt2, 'FontSize', 9);

%% Modulation index
frB = abs(FR_baseline_all);
frE = abs(FR_eopn_all);

% Modulation index (light − control) / (light + control)
modIdx = (frE - frB) ./ (frE + frB);   % attention-style index [web:64][web:65]

% Optional: handle 0/0 or tiny denominators
modIdx(abs(frE + frB) < 1e-6) = NaN;
[p,h,stats] = signrank(modIdx);    

mi_mean = nanmean(modIdx);
mi_median = nanmedian(modIdx);

figure; hold on

h = histogram(modIdx, 'NumBins', 20, ...
                       'FaceColor', [0.1 0.1 0.1], ...
                       'EdgeColor', 'none');            % [web:88]

xline(0,'--','Color',[0.6 0.6 0.6],'LineWidth',1.5);    % zero line [web:68]

xlabel('Change in modulation index (control-light)');
ylabel('Neurons');
title('eOPN3');

% Place stats text near top-right of axes [web:78][web:83]
yl = ylim;
xl = xlim;
txt = sprintf('median = %.2f\np = %.3g (Wilcoxon)', mi_median, p);
text(xl(1)+0.55*range(xl), yl(1)+0.9*range(yl), txt, ...
     'FontSize', 9, 'HorizontalAlignment','left');
xlim([-1 1])
set(gca,'Box','off','TickDir','out','FontSize',9);axis square

%%

%%
% index of unit with largest absolute modulation
[~,bestIdx] = max(abs(modIdx));   % same ordering as FR_*_all
modIdx2 = modIdx(end-length(baselineSpks):end);
[~,bestIdx] = max(abs(modIdx2));   % same ordering as FR_*_all
bestIdx = 21;
% If you want to inspect a specific session, adapt indices accordingly.
B_unit = baselineSpks{bestIdx};   % [nTrials x nTime]
E_unit = eopnSpks{bestIdx};       % [nTrials x nTime]

t = (1:size(B_unit,2));           % ms bins relative to trial start
stimOn  = win1;                   % ms index when light turns on
stimOff = win2;                   % ms index when light turns off
figure;

% ----- Top subplot: raster -----
ax1 = subplot(2,1,1); hold on

% baseline trials first (rows 1..nB)
[nB,~] = size(B_unit);
[nE,~] = size(E_unit);

% baseline spikes: black dots
[trialIdx_B, timeIdx_B] = find(B_unit);   % [row,col] of 1s
scatter(t(timeIdx_B), trialIdx_B, 8, 'k', 'filled');           % [web:6]

% eOPN spikes: magenta dots, stacked after baseline
[trialIdx_E, timeIdx_E] = find(E_unit);
scatter(t(timeIdx_E), trialIdx_E + nB, 8, [0.6 0 0.6], 'filled');

% light window as vertical band or lines
xline(stimOn,'--','Color',[0.5 0.5 0.5]);                      % [web:68]
xline(stimOff,'--','Color',[0.5 0.5 0.5]);
ylabel('Trials');
title('Top modulated unit');

ylim([0 nB+nE+1]);
ylim([0 90])
set(gca,'YDir','reverse','Box','off','TickDir','out');         % raster style [web:23]
% bin for PSTH (could be same or smaller than analysis bin)
psthBin = 20;                        % ms
edges = 1:psthBin:(size(B_unit,2)+1);
binCenters = edges(1:end-1) + psthBin/2;

% baseline PSTH
cntB = zeros(nB, numel(edges)-1);
for b = 1:numel(edges)-1
    idx = edges(b):edges(b+1)-1;
    cntB(:,b) = sum(B_unit(:,idx),2);
end
rateB = mean(cntB,1) * (1000/psthBin);   % Hz [web:19]

% eOPN PSTH
cntE = zeros(nE, numel(edges)-1);
for b = 1:numel(edges)-1
    idx = edges(b):edges(b+1)-1;
    cntE(:,b) = sum(E_unit(:,idx),2);
end
rateE = mean(cntE,1) * (1000/psthBin);   % Hz [web:19]

% ----- Bottom subplot: PSTH -----
ax2 = subplot(2,1,2); hold on
plot(binCenters, rateB, 'Color',[0.4 0.4 0.4], 'LineWidth',1.5);       % baseline (gray)
plot(binCenters, rateE, 'Color',[0.8 0 0.6],   'LineWidth',1.5);       % eOPN (magenta)

xline(stimOn,'--','Color',[0.5 0.5 0.5]);                              % light window [web:68]
xline(stimOff,'--','Color',[0.5 0.5 0.5]);

xlabel('Time (ms)');
ylabel('Firing rate (spikes/s)');
set(gca,'Box','off','TickDir','out');

% share x-axis between raster and PSTH
linkaxes([ax1,ax2],'x');                                              % [web:98]


function output = make_nice_mean_raster(spmat,smooth_window,showplot)
%*********** spmat1 and spmat2 are spike matrices of two conditions you wish to compare
%*********** smooth_window ... gaussian smoothing in millisecs
numconds = size(spmat,2);
if (numconds==2)
    colo = [[1,0,0];[0,0,1]];
else
    colo = jet(numconds);
end
for k = 1:numconds
    spud = spmat{k};
    numtrials = size(spud,1);
    smorate = gauss_smooth(sum( spud(1:numtrials,:))/....
        numtrials,smooth_window)*1000;
    if showplot
        plot(smorate,'k'); hold on;
        %                 set(H,'Color',colo(k,:));
    end
    output(k,:) = smorate;

end
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

function [normSpikeRate,idx,idxc] = spknorm(temp)
[nanIdx,~,~] = find(~isnan(temp));
nanIdx = unique(nanIdx);
normSpikeRate = zscore(temp(nanIdx,:),0,2);
idx = zeros(size(normSpikeRate,1),1);
for n = 1:length(idx)
    [~,idx(n)] = max(normSpikeRate(n,:));
end
[~,idxc] = sort(idx);
end
function spk_opto = adOpto(tempSpk,optoIdx)
spk_mod = tempSpk;        % copy to modify


optoFreq    = 20;         % Hz
binSize_ms  = 1;          % your current spike bin size
optoStart   = 1580;       % first bin to consider (ms)
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

function plotSpkSeq(normSpikeRate, idxc)
idx = zeros(size(normSpikeRate,1),1);
for n = 1:length(idx)
    [~,idx(n)] = max(normSpikeRate(n,:));
end
if ~exist('idxc','var')
    [~,idxc] = sort(idx);
end
path = idx(idxc);

% % Mean FR for each neuron (across time)
% meanFR = max(spikeRate,[], 2);          % size: [neurons x 1]
% meanFR_sorted = meanFR(idxc);             % reorder by idxc

% Create two axes: heatmap and mean FR
% figure;
% ax1 = subplot(1,2,1);                      % left: heatmap
imagesc(-0.5*1000:1.5*1000, ...
    1:size(normSpikeRate,1), ...
    normSpikeRate(idxc,:));
hold on;
plot((path)-0.5*1000, 1:size(normSpikeRate,1), 'r', 'LineWidth', 1);
xlabel('Time (ms)');
ylabel('Neuron (sorted)');
caxis([0.0 2])
% xlim([-500 1500])

% ax2 = subplot(1,2,2);                      % right: mean FR
% barh(1:size(meanFR_sorted,1),meanFR_sorted,'k');
% set(ax2, 'YDir', 'reverse');               % match imagesc orientation
% ylim([0.5 size(spikeRate,1)+0.5]);
% xlabel('Mean FR');
% yticklabels([]);                           % hide duplicate y labels
% linkaxes([ax1 ax2],'y');                   % keep neuron order aligned
% axis off
end