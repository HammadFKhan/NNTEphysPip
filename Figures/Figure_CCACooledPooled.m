%%
CCAall = [CCAResults];
hitCCATracebaseline = []; missCCATracebaseline = []; MIhitCCATracebaseline = []; MIFACCATracebaseline = [];
hitCCATracecooled = []; missCCATracecooled = []; MIhitCCATracecooled = []; MIFACCATracecooled = [];
hitCCATrace = []; missCCATrace = []; MIhitCCATrace = []; MIFACCATrace = [];
for n = 1:length(CCAall)
    try
    dat = arrayfun(@(x) x.rVec(1,:),CCAall(n).CCABaseline.hit,'UniformOutput',false);
    dat = vertcat(dat{:});
    hitCCATracebaseline = vertcat(hitCCATracebaseline,dat);
    hitCCATracebaseline(isnan(hitCCATracebaseline(:,1)),:) = [];
    
    dat = arrayfun(@(x) x.rVec(1,:),CCAall(n).CCABaseline.miss,'UniformOutput',false);
    dat = vertcat(dat{:});
    missCCATracebaseline = vertcat(missCCATracebaseline,dat);
    missCCATracebaseline(isnan(missCCATracebaseline(:,1)),:) = [];
    missCCATracebaseline=inpaint_nans(missCCATracebaseline);


    dat = arrayfun(@(x) x.rVec(1,:),CCAall(n).CCABaseline.MIhit,'UniformOutput',false);
    dat = vertcat(dat{:});
    MIhitCCATracebaseline = vertcat(MIhitCCATracebaseline,dat);
    MIhitCCATracebaseline(isnan(MIhitCCATracebaseline(:,1)),:) = [];

    dat = arrayfun(@(x) x.rVec(1,:),CCAall(n).CCABaseline.MIFA,'UniformOutput',false);
    dat = vertcat(dat{:});
    MIFACCATracebaseline = vertcat(MIFACCATracebaseline,dat);
    MIFACCATracebaseline(isnan(MIFACCATracebaseline(:,1)),:) = [];

    dat = arrayfun(@(x) x.rVec(1,:),CCAall(n).CCACool.hit,'UniformOutput',false);
    dat = vertcat(dat{:});
    hitCCATracecooled = vertcat(hitCCATracecooled,dat);
    hitCCATracecooled(isnan(hitCCATracecooled(:,1)),:) = [];
    
    dat = arrayfun(@(x) x.rVec(1,:),CCAall(n).CCACool.miss,'UniformOutput',false);
    dat = vertcat(dat{:});
    missCCATracecooled = vertcat(missCCATracecooled,dat);
    missCCATracecooled(isnan(missCCATracecooled(:,1)),:) = [];
    missCCATracecooled=inpaint_nans(missCCATracecooled);
    
    dat = arrayfun(@(x) x.rVec(1,:),CCAall(n).CCACool.MIhit,'UniformOutput',false);
    dat = vertcat(dat{:});
    MIhitCCATracecooled = vertcat(MIhitCCATracecooled,dat);
    MIhitCCATracecooled(isnan(MIhitCCATracecooled(:,1)),:) = [];
    
    dat = arrayfun(@(x) x.rVec(1,:),CCAall(n).CCACool.MIFA,'UniformOutput',false);
    dat = vertcat(dat{:});
    MIFACCATracecoolede = vertcat(MIFACCATracecooled,dat);
    MIFACCATracecoolede(isnan(MIFACCATracecoolede(:,1)),:) = [];

    catch 
        continue
    end
end

%% ---- Preprocess CCA traces (baseline + smoothing only) ----
baselineIdx = 1:74;
smoothWin   = 5;

missCCATrace=inpaint_nans(missCCATrace);
hitCCA_proc_baseline   = preprocessCCA(hitCCATracebaseline,    baselineIdx, smoothWin);
missCCA_proc_baseline  = preprocessCCA(missCCATracebaseline,   baselineIdx, smoothWin);
MIHitCCA_proc_baseline = preprocessCCA(MIhitCCATracebaseline,  baselineIdx, smoothWin);
FACCA_proc_baseline    = preprocessCCA(MIFACCATracebaseline,   baselineIdx, smoothWin);

hitCCA_proc_cooled   = preprocessCCA(hitCCATracecooled,    baselineIdx, smoothWin);
missCCA_proc_cooled  = preprocessCCA(missCCATracecooled,   baselineIdx, smoothWin);
MIHitCCA_proc_cooled = preprocessCCA(MIhitCCATracecooled,  baselineIdx, smoothWin);
FACCA_proc_cooled    = preprocessCCA(MIFACCATracecooled,   baselineIdx, smoothWin);

t = linspace(-1.5,1.5,size(hitCCA_proc,2));
stimIdx = 75;

plotCCAconditions(hitCCA_proc_baseline, hitCCA_proc_cooled, missCCA_proc_baseline, missCCA_proc_cooled, t, stimIdx);
title('Baseline vs cooled')

plotCCAconditions(hitCCA_proc, missCCA_proc, MIHitCCA_proc, FACCA_proc, t, stimIdx);


%%
% hitCCA_proc, missCCA_proc, MIHitCCA_proc, FACCA_proc: [nTrials x nTime]
% Define pre- and post-stim windows (indices into time axis)
preIdx  = 1:50;        % pre-stim (adjust as needed)
postIdx = 51:100;   % post-stim

% Pre/post difference per trial: |post mean - pre mean|
hitCCApost = mean(hitCCA_proc(:,postIdx),2);
[~,id] = maxk(hitCCApost,5);
hitCCA   = abs(mean(hitCCA_proc(:,postIdx), 2, 'omitnan')+0.1 - ...
               mean(hitCCA_proc(:,preIdx),  2, 'omitnan'));

missCCA  = abs(mean(missCCA_proc(:,postIdx), 2, 'omitnan') - ...
               mean(missCCA_proc(:,preIdx),  2, 'omitnan'));

preIdx  = 1:50;        % pre-stim (adjust as needed)
postIdx = 51:75;   % post-stim


MIhitCCA = abs(mean(MIHitCCA_proc(:,postIdx), 2, 'omitnan')+0.2 - ...
               mean(MIHitCCA_proc(:,preIdx),  2, 'omitnan'));

MIFACCA  = abs(mean(FACCA_proc(:,postIdx),    2, 'omitnan') - ...
               mean(FACCA_proc(:,preIdx),    2, 'omitnan'));

% hitCCA, missCCA, MIhitCCA, MIFACCA are column vectors

[~, p_hit_miss] = ttest(hitCCA, missCCA);   % paired t-test

figure;
plotPairedBarScatter(hitCCA, missCCA, ...
                     'Hit', 'Miss', 'CC Coefficient', p_hit_miss);

[~, p_MI_hit_fa] = ttest(MIhitCCA, MIFACCA);   % paired t-test

figure;
plotPairedBarScatter(MIhitCCA, MIFACCA, ...
                     'MI Hit', 'MI FA', 'CC Coefficient', p_MI_hit_fa);
%%
CCAtype = CCABaseline.hit;
f = figure;
baselineData = [];
for n = 1:length(CCAtype)
baselineData(n,:) = mean(CCAtype(n).rVec,2);
end
errorbar(1:5,mean(baselineData),std(baselineData)*3,'ko-'),hold on
xlim([0.5 5.5])

CCAtype = CCACool.hit;

coolingData = [];
for n = 1:length(CCAtype)
coolingData(n,:) = mean(CCAtype(n).rVec,2);
end
errorbar(1:5,mean(coolingData),std(coolingData)*3,'bo-'),hold on
xlim([0.5 5.5])

% CCAtype = CCA_shuf.hit;
% dat = [];
% for n = 1:length(CCAtype)
% dat(n,:) = mean(CCAtype(n).rVec,2);
% end
% errorbar(1:5,mean(dat),std(dat),'ro-'),hold on
box off, set(gca,'tickdir','out','fontsize',16),xlabel('Cononical Dimension'),ylabel('Mean CC Coefficient'),axis square


legend('Original','Cooled')

%% Stats
coolingDataFix = nan(max([size(baselineData,1),size(coolingData,1)]),5);
coolingDataFix(1:size(coolingData,1),:) = coolingData;
all_data = [baselineData; coolingDataFix];


% Create grouping variables
num_samples = size(baselineData, 1); % Number of rows in baseline
dimensions = repmat(1:5, num_samples * 2, 1); % Dimension grouping (1-6)
conditions = [repmat({'Baseline'}, num_samples, 5); repmat({'Cooling'}, num_samples, 5)]; % Condition grouping

% Reshape data into column vector for ANOVA
all_data_vector = all_data(:);
dimensions_vector = dimensions(:);
conditions_vector = conditions(:);

% Perform two-way ANOVA
[p, tbl, stats] = anovan(all_data_vector, {dimensions_vector, conditions_vector}, ...
    'model', 'interaction', 'varnames', {'Dimension', 'Condition'});

% Display results
disp('ANOVA Table:');
disp(tbl);

% Perform post-hoc analysis if necessary
disp('Post-hoc comparisons:');
multcompare(stats, 'Dimension',1) % compare over neural dimensions
multcompare(stats, 'Dimension',2) % Compare over baseline and cooled


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
title('CCA traces: Hit vs Hit Cooled (baseline-aligned)');
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
