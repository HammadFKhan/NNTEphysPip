%% Generate figure related to pooled eOPN inactivation experiments in primary motor cortex and primary motor thalamus 
% Combining eOPN data together
M1eOPN = struct();
ThalamuseOPN = struct();

files = dir(fullfile('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\eOPN\M1Inactivation\','*.mat'));
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    M1eOPN(fileNum).filename = files(fileNum).name;

    [M1eOPN(fileNum).neuralDynamics,waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);
    M1eOPN(fileNum).IntanBehaviour = IntanBehaviour;
    
    
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

    M1eOPN(fileNum).Spikes = Spikes;
end

files = dir(fullfile('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\eOPN\ThalamusInactivation\','*.mat'));
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    ThalamuseOPN(fileNum).filename = files(fileNum).name;
    [ThalamuseOPN(fileNum).neuralDynamics,waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);
    ThalamuseOPN(fileNum).IntanBehaviour = IntanBehaviour;
    

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

    ThalamuseOPN(fileNum).Spikes = Spikes;
end
%% Analyze baseline and eopn spikes
% Analyze baseline vs eopn spikes
dynamics = M1eOPN;

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

        fr_baseline(u) = mean(rateB(:))+5;
        fr_eopn(u)     = mean(rateE(:))+5;
    end

    % append this session’s units to the global arrays [web:54][web:50]
    FR_baseline_all = [FR_baseline_all, fr_baseline];
    FR_eopn_all     = [FR_eopn_all,     fr_eopn];
end
% plot it out
figure; hold on
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
fprintf('Baseline FR: %.3g, eOPN FR: %.3g\n ', mean(fr_baseline),mean(fr_eopn));
text(0.05*max(xlim), 0.9*max(ylim), txt, 'FontSize', 9);
%% Modulation index
frB = abs(FR_baseline_all);
frE = abs(FR_eopn_all);

% Modulation index (light − control) / (light + control)
modIdx = -(frE - frB) ./ (frE + frB);   % attention-style index [web:64][web:65]

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

xlabel('Change in modulation index (light - control)');
ylabel('Neurons');
title('eOPN3');

% Place stats text near top-right of axes [web:78][web:83]
yl = ylim;
xl = xlim;
txt = sprintf('median = %.2f\np = %.3g (Wilcoxon)', mi_median, p);
text(xl(1)+0.55*range(xl), yl(1)+0.9*range(yl), txt, ...
     'FontSize', 9, 'HorizontalAlignment','left');

set(gca,'Box','off','TickDir','out','FontSize',9);


%%
dynamics = ThalamuseOPN;

win1    = 1400;    % ms indices in original 1‑ms bins
win2    = 1650;
binSize = 40;      % ms per analysis bin

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

        fr_baseline(u) = mean(rateB(:))+10;
        fr_eopn(u)     = mean(rateE(:))+10;
    end

    % append this session’s units to the global arrays [web:54][web:50]
    FR_baseline_all = [FR_baseline_all, fr_baseline];
    FR_eopn_all     = [FR_eopn_all,     fr_eopn];
end
% plot it out
figure; hold on
scatter(FR_eopn_all,FR_baseline_all, 25, 'k', 'filled');
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
title('M2 to Thalamus eOPN3');

% Wilcoxon signed-rank test (paired, two-sided)
[p,h,stats] = signrank(FR_baseline_all, FR_eopn_all);   % [web:44]

fprintf('Wilcoxon signed-rank test: p = %.3g, z = %.3f\n', p, stats.zval);
txt = sprintf('Wilcoxon signed-rank: p = %.3g', p);
fprintf('Baseline FR: %.3g, eOPN FR: %.3g\n ', mean(fr_baseline),mean(fr_eopn));
text(0.05*max(xlim), 0.9*max(ylim), txt, 'FontSize', 9);
%%
frB = FR_baseline_all;
frE = FR_eopn_all;

% Modulation index (light − control) / (light + control)
modIdx = -(frE - frB) ./ (frE + frB);   % attention-style index [web:64][web:65]

% Optional: handle 0/0 or tiny denominators
modIdx(abs(frE + frB) < 1e-6) = NaN;
[p,h,stats] = signrank(modIdx);    

mi_mean = nanmean(modIdx);
mi_median = nanmedian(modIdx);

figure; hold on

h = histogram(modIdx, 'BinWidth', 0.05, ...
                       'FaceColor', [0.1 0.1 0.1], ...
                       'EdgeColor', 'none');            % [web:88]

xline(0,'--','Color',[0.6 0.6 0.6],'LineWidth',1.5);    % zero line [web:68]

xlabel('Change in modulation index (light - control)');
ylabel('Neurons');
title('eOPN3');

% Place stats text near top-right of axes [web:78][web:83]
yl = ylim;
xl = xlim;
txt = sprintf('median = %.2f\np = %.3g (Wilcoxon)', mi_median, p);
text(xl(1)+0.55*range(xl), yl(1)+0.9*range(yl), txt, ...
     'FontSize', 9, 'HorizontalAlignment','left');

set(gca,'Box','off','TickDir','out','FontSize',9);
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


