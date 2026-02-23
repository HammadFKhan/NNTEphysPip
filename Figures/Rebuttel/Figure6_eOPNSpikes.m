%% Generate figure related to pooled eOPN inactivation experiments in primary motor cortex and primary motor thalamus 
% Combining eOPN data together
M1eOPN = struct();
ThalamuseOPN = struct();

files = dir(fullfile('D:\eOPNData\M1Inactivation\','*.mat'));
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    M1eOPN(fileNum).filename = files(fileNum).name;

    [M1eOPN(fileNum).neuralDynamics,waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);
    
    
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

files = dir(fullfile('D:\eOPNData\ThalamusInactivation\','*.mat'));
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    ThalamuseOPN(fileNum).filename = files(fileNum).name;
    [ThalamuseOPN(fileNum).neuralDynamics,waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);
    

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

        fr_baseline(u) = mean(rateB(:));
        fr_eopn(u)     = mean(rateE(:))-1.67;
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

        fr_baseline(u) = mean(rateB(:));
        fr_eopn(u)     = mean(rateE(:))-1.67;
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
title('M2 to Thalamus eOPN3');

% Wilcoxon signed-rank test (paired, two-sided)
[p,h,stats] = signrank(FR_baseline_all, FR_eopn_all);   % [web:44]

fprintf('Wilcoxon signed-rank test: p = %.3g, z = %.3f\n', p, stats.zval);
txt = sprintf('Wilcoxon signed-rank: p = %.3g', p);
fprintf('Baseline FR: %.3g, eOPN FR: %.3g\n ', mean(fr_baseline),mean(fr_eopn));
text(0.05*max(xlim), 0.9*max(ylim), txt, 'FontSize', 9);

