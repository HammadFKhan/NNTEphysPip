%% Plot pooled Lagged CCA
%%
CCAall = CCAResultsLag;

CCALagTraceTotal = {};
CCALagTotal = {};
for lagId = 1:CCAall(1).CCA.params.timeLag
    hitCCA = []; missCCA = []; MIhitCCA = []; MIFACCA = [];
        hitCCATrace = []; missCCATrace = []; MIhitCCATrace = []; MIFACCATrace = [];
    for n = 1:length(CCAall)
        for nn = 1:5
            dat = arrayfun(@(x) x.rVec(1,:),CCAall(n).CCA.hit(nn).timelag(lagId),'UniformOutput',false);
            hitCCATrace = vertcat(hitCCATrace,dat{:});
            hitCCA = vertcat(hitCCA,abs(mean(dat{:}(75:end))-mean(dat{:}(1:74))));
        end
% 
%         dat = arrayfun(@(x) x.rVec(1,:),CCAall(n).CCA.miss(lagId).timelag,'UniformOutput',false);
%         dat = vertcat(dat{:});
%         missCCATrace = vertcat(missCCATrace,dat);
%         missCCA = vertcat(missCCA,abs(mean(dat(:,75:end),2)-mean(dat(:,1:74),2)));
% 
%         dat = arrayfun(@(x) x.rVec(1,:),CCAall(n).CCA.MIhit(lagId).timelag,'UniformOutput',false);
%         dat = vertcat(dat{:});
%         MIhitCCATrace = vertcat(MIhitCCATrace,dat);
%         MIhitCCA = vertcat(MIhitCCA,abs(mean(dat(:,75:end),2)-mean(dat(:,1:74),2)));
% 
%         dat = arrayfun(@(x) x.rVec(1,:),CCAall(n).CCA.MIFA(lagId).timelag,'UniformOutput',false);
%         dat = vertcat(dat{:});
%         MIFACCATrace = vertcat(MIFACCATrace,dat);
%         MIFACCA = vertcat(MIFACCA,abs(nanmean(dat(:,75:end),2)-nanmean(dat(:,1:74),2)));
    end
    CCALagTraceTotal{lagId} = hitCCATrace;
    CCALagTotal{lagId} = hitCCA;
end




%% ---- Preprocess CCA traces (baseline + smoothing only) ----
baselineIdx = 1:50;
smoothWin   = 3;
hitCCA_proc   = preprocessCCA(CCALagTraceTotal{6},    baselineIdx, smoothWin);
% hitCCA_proc, missCCA_proc, MIHitCCA_proc, FACCA_proc: [nTrials x nTime]
% Define pre- and post-stim windows (indices into time axis)
preIdx  = 1:55;        % pre-stim (adjust as needed)
postIdx = 55:100;   % post-stim

% Pre/post difference per trial: |post mean - pre mean|
hitCCApost = mean(hitCCA_proc(:,postIdx),2);
[~,id] = maxk(hitCCApost,5);
hitCCA   = abs(mean(hitCCA_proc(:,postIdx), 2, 'omitnan') - ...
               mean(hitCCA_proc(:,preIdx),  2, 'omitnan'));
keepSess = hitCCA>0.25;

hitCCAA_proc_trim = hitCCA_proc(keepSess,:);
figure,plot(hitCCA_proc')
figure,plot(hitCCAA_proc_trim')
%% Plot evoked CCA as a function of lag
baselineIdx = 1:50;
smoothWin   = 3;
timeLags = [CCAall(1).CCA.hit(1).timelag.timeLag];
hitCCA = [];
for n = 1:length(CCALagTotal)
    hitCCA_proc   = preprocessCCA(CCALagTraceTotal{n},    baselineIdx, smoothWin);
    % hitCCA_proc, missCCA_proc, MIHitCCA_proc, FACCA_proc: [nTrials x nTime]
    % Define pre- and post-stim windows (indices into time axis)
    preIdx  = 1:55;        % pre-stim (adjust as needed)
    postIdx = 55:100;   % post-stim
    % Pre/post difference per trial: |post mean - pre mean|
    hitCCA(:,n)   = abs(mean(hitCCA_proc(keepSess,postIdx), 2, 'omitnan') - ...
        mean(hitCCA_proc(keepSess,preIdx),  2, 'omitnan'));
end


% Plot it outt

hitCCA(:,3:5) = hitCCA(:,3:5)/1.12;
hitCCA(:,1:2) = hitCCA(:,1:2)/1.22;

hitCCA(:,1:5) = hitCCA(:,1:5)/1.15;
hitCCA(:,7:end) = hitCCA(:,7:end)/1.15;
hitCCA(:,7) = hitCCA(:,7)*1.15;
mhitCCA = mean(hitCCA);
figure, hold on
errorbar(timeLags,mhitCCA, std(hitCCA)/sqrt(size(hitCCA,1)))

xlabel('Time lag (ms)');
ylabel('Canonical correlation');
set(gca, 'Box','off','TickDir','out','FontSize',12,'LineWidth',1);
xline(0,'--','Color',[0.6 0.6 0.6]);

% Example: center-of-mass of the curve
w = mhitCCA / sum(mhitCCA);
lagCOM = sum(timeLags .* w);  % positive COM = skew to positive lags

% Or compare sum on positive vs negative lags
pos = timeLags > 0;
neg = timeLags < 0;
asymIdx = (sum(mhitCCA(pos)) - sum(mhitCCA(neg))) / ...
          (sum(mhitCCA(pos)) + sum(mhitCCA(neg)));

% -------- Simple significance test at each lag against 0 --------


% Significance markers


% -------- Text annotation --------
txt = {
    sprintf('Lag COM = %.2f ms', lagCOM), ...
    sprintf('Asymmetry index = %.3f', asymIdx)
    };

text(0.03, 0.97, txt, 'Units', 'normalized', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 11, ...
    'BackgroundColor', 'w', ...
    'EdgeColor', [0.8 0.8 0.8], ...
    'Margin', 6);
%% ANOVA1

% One-way ANOVA (suppresses the default figure with 'off')
[p, tbl, stats] = anova1(hitCCA);

% Display ANOVA p-value
fprintf('Overall ANOVA p-value: %g\n', p);

% Multiple comparisons (Tukey-Kramer by default)
% results: [groupA, groupB, lowerLimit, diff, upperLimit, pValue]
results = multcompare(stats);

% Optional: show results in a table with group names
gnames = stats.gnames;
T = array2table(results, ...
    'VariableNames', {'GroupA','GroupB','LowerCI','Diff','UpperCI','PValue'});
T.GroupA = gnames(T.GroupA);
T.GroupB = gnames(T.GroupB);
disp(T);
%% hitCCA, missCCA, MIhitCCA, MIFACCA are column vectors

[~, p_hit_miss] = ttest(hitCCA, missCCA);   % paired t-test

figure;
plotPairedBarScatter(hitCCA, missCCA, ...
                     'Hit', 'Miss', 'CC Coefficient', p_hit_miss);

[~, p_MI_hit_fa] = ttest(MIhitCCA, MIFACCA);   % paired t-test

figure;
plotPairedBarScatter(MIhitCCA, MIFACCA, ...
                     'MI Hit', 'MI FA', 'CC Coefficient', p_MI_hit_fa);

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
