%% Generate figure related to speed differences based on matched firing rates
% Combining eOPN data together
M1eOPN = struct();
ThalamuseOPN = struct();

files = dir(fullfile('D:\eOPNData\M1Inactivation\','*.mat'));
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    % Calculated mean matched trajectory speed across conditions
    [M1eOPN(fileNum).matchedData] = meanRateSpeed(Spikes,IntanBehaviour);
    M1eOPN(fileNum).filename = files(fileNum).name;
end

files = dir(fullfile('D:\eOPNData\ThalamusInactivation\','*.mat'));
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    [ThalamuseOPN(fileNum).matchedData] = meanRateSpeed(Spikes,IntanBehaviour);
    ThalamuseOPN(fileNum).filename = files(fileNum).name;
end
%% Plot matched trajectory and spike dynamics
% Decompose example session
fileNum = 2;
valid = M1eOPN(fileNum).matchedData.valid;
frValid = M1eOPN(fileNum).matchedData.frValid;
condValid = M1eOPN(fileNum).matchedData.condValid;
spdValid = M1eOPN(fileNum).matchedData.spdValid;
edges = linspace(prctile(frValid, 5), prctile(frValid, 95), 20);
[~, ~, binIdx] = histcounts(frValid, edges);

nBins = max(binIdx);
meanSpeed_control = nan(nBins, 1);
meanSpeed_manip   = nan(nBins, 1);

for b = 1:nBins
    inBin = binIdx == b;
    inControl = inBin & (condValid == 0);
    inManip   = inBin & (condValid == 1);
    meanSpeed_control(b) = mean(spdValid(inControl));
    meanSpeed_manip(b)   = nanmean(spdValid(inManip));
end

% Min-max normalize firing rate and speed to [0,1]
frMin = min(frValid);
frMax = max(frValid);
frNorm = (frValid - frMin) / (frMax - frMin);

spdMin = min(spdValid);
spdMax = max(spdValid);
spdNorm = (spdValid - spdMin) / (spdMax - spdMin);

% Bin on raw mean firing rate
%     edges = linspace(prctile(frValid, 5), prctile(frValid, 95), 20);  % 20 bins
edges = linspace(0, 1, 20);  % 20 bins
[~, ~, binIdx] = histcounts(frNorm, edges);

nBins = max(binIdx);
meanSpeed_control = nan(nBins, 1);
meanSpeed_manip   = nan(nBins, 1);

for b = 1:nBins
    inBin   = binIdx == b;
    inCtrl  = inBin & (condValid == 0);
    inManip = inBin & (condValid == 1);

    meanSpeed_control(b) = nanmean(spdNorm(inCtrl));
    meanSpeed_manip(b)   = mean(spdNorm(inManip));
end

binCenters = 0.5 * (edges(1:end-1) + edges(2:end));

% ---- Plot ----
figure; hold on;
plot(binCenters, meanSpeed_control, '-o');
plot(binCenters, meanSpeed_manip,   '-o');
xlabel('Normalized mean firing rate');
ylabel('Normalized trajectory speed');
legend({'Control','Manip'}, 'Location','best');
title('Speed vs mean firing (normalized, matched bins)');

grayLight = [0.7 0.7 0.7];
grayDark  = [0.3 0.3 0.3];
orange    = [0.85 0.45 0.1];  % adjust to match your scheme

% Assume you already have: binCenters, meanSpeed_control, meanSpeed_manip

validCtrl  = ~isnan(meanSpeed_control);
validManip = ~isnan(meanSpeed_manip);

x_ctrl  = binCenters(validCtrl)';
y_ctrl  = meanSpeed_control(validCtrl);
x_manip = binCenters(validManip)';
y_manip = meanSpeed_manip(validManip);

% Fit straight lines to binned data
p_ctrl  = polyfit(x_ctrl,  y_ctrl,  1);
p_manip = polyfit(x_manip, y_manip, 1);

x_fit = linspace(min(binCenters), max(binCenters), 100);
y_fit_ctrl  = polyval(p_ctrl,  x_fit);
y_fit_manip = polyval(p_manip, x_fit);

figure; hold on;

% Binned points
plot(x_ctrl,  y_ctrl,  'o', 'MarkerFaceColor', grayLight, ...
    'MarkerEdgeColor', grayDark);
plot(x_manip, y_manip, 'o', 'MarkerFaceColor', orange, ...
    'MarkerEdgeColor', orange);

% Fitted lines
plot(x_fit, y_fit_ctrl,  '-', 'Color', grayDark, 'LineWidth', 2);
plot(x_fit, y_fit_manip, '-', 'Color', orange,   'LineWidth', 2);

xlabel('Mean firing rate');
ylabel('Trajectory speed');
legend({'Baseline bins','eOPN bins','Baseline fit','eOPN fit'}, ...
    'Location','best');
title('Speed vs mean firing (binned data with linear fits)');
box off;
%% Loop over session
optoType = M1eOPN;
nSessions = length(optoType);
nBins     = 20;
edges     = linspace(0,1,nBins+1);
binCenters = 0.5*(edges(1:end-1)+edges(2:end));

% Store per-session binned curves
all_ctrl  = nan(nSessions, nBins);
all_manip = nan(nSessions, nBins);

for fileNum = 1:nSessions
    valid    = optoType(fileNum).matchedData.valid;
    frValid  = optoType(fileNum).matchedData.frValid(valid);
    spdValid = optoType(fileNum).matchedData.spdValid(valid);
    condValid = optoType(fileNum).matchedData.condValid(valid);

    % Min-max normalize within session
    frMin = min(frValid);
    frMax = max(frValid);
    frNorm = (frValid - frMin) / (frMax - frMin);

    spdMin = min(spdValid);
    spdMax = max(spdValid);
    spdNorm = (spdValid - spdMin) / (spdMax - spdMin);

    % Bin on normalized firing rate
    [~, ~, binIdx] = histcounts(frNorm, edges);
    nBinsSess = max(binIdx);  % could be < nBins if outer bins empty

    meanSpeed_control = nan(1, nBins);
    meanSpeed_manip   = nan(1, nBins);

    for b = 1:nBinsSess
        inBin   = binIdx == b;
        inCtrl  = inBin & (condValid == 0);
        inManip = inBin & (condValid == 1);

        if any(inCtrl)
            meanSpeed_control(b) = mean(spdNorm(inCtrl));
        end
        if any(inManip)
            meanSpeed_manip(b) = mean(spdNorm(inManip));
        end
    end

    all_ctrl(fileNum, :)  = meanSpeed_control;
    all_manip(fileNum, :) = meanSpeed_manip;
end

% Session-averaged curves and SEM
mean_ctrl = nanmean(all_ctrl, 1);
mean_manip = nanmean(all_manip, 1);

sem_ctrl  = nanstd(all_ctrl, [], 1) ./ sqrt(sum(~isnan(all_ctrl),1));
sem_manip = nanstd(all_manip, [], 1) ./ sqrt(sum(~isnan(all_manip),1));

grayLight = [0.7 0.7 0.7];
grayDark  = [0.3 0.3 0.3];
orange    = [0.85 0.45 0.1];
% Plot session-averaged binned curves with SEM
figure; hold on;

% Control
errorbar(binCenters, mean_ctrl, sem_ctrl, '-o', ...
    'Color', grayDark, 'MarkerFaceColor', grayLight, ...
    'MarkerEdgeColor', grayDark);

% Manip
errorbar(binCenters, mean_manip, sem_manip, '-o', ...
    'Color', orange, 'MarkerFaceColor', orange, ...
    'MarkerEdgeColor', orange);

xlabel('Normalized mean firing rate');
ylabel('Normalized trajectory speed');
legend({'Baseline','eOPN'}, 'Location','best');
title('Speed vs mean firing across sessions (normalized, binned)');
set(gca,'tickdir','out','fontsize',10)
box off;
axis square

% ---------- SESSION-AVERAGED BINS ----------
mean_ctrl   = nanmean(all_ctrl, 1);
mean_manip  = nanmean(all_manip, 1);

validCtrl   = ~isnan(mean_ctrl);
validManip  = ~isnan(mean_manip);

x_ctrl  = binCenters(validCtrl);
y_ctrl  = mean_ctrl(validCtrl);
x_manip = binCenters(validManip);
y_manip = mean_manip(validManip);

% ---------- FIT + CI USING YOUR CODE ----------
grayLight = [0.7 0.7 0.7];
grayDark  = [0.3 0.3 0.3];
orange    = [0.85 0.45 0.1];

% Control: fit linear model to binned points
tbl_ctrl = table(x_ctrl', y_ctrl', 'VariableNames', {'FR','Speed'});
mdl_ctrl = fitlm(tbl_ctrl, 'Speed ~ FR');

% Manip: fit linear model
tbl_manip = table(x_manip', y_manip', 'VariableNames', {'FR','Speed'});
mdl_manip = fitlm(tbl_manip, 'Speed ~ FR');

% Grid for plotting fits and CIs
x_fit = linspace(min(binCenters), max(binCenters), 100)';

[y_ctrl_fit,  y_ctrl_ci]  = predict(mdl_ctrl,  table(x_fit,'VariableNames',{'FR'}), 'Prediction','curve');
[y_manip_fit, y_manip_ci] = predict(mdl_manip, table(x_fit,'VariableNames',{'FR'}), 'Prediction','curve');

figure; hold on;

% Binned session-averaged points
plot(x_ctrl,  y_ctrl,  'o', 'MarkerFaceColor', grayLight, ...
    'MarkerEdgeColor', grayDark);
plot(x_manip, y_manip, 'o', 'MarkerFaceColor', orange, ...
    'MarkerEdgeColor', orange);
axis square;

% Confidence bands (patches)
ctrl_fillX = [x_fit; flipud(x_fit)];
ctrl_fillY = [y_ctrl_ci(:,1); flipud(y_ctrl_ci(:,2))];
patch(ctrl_fillX, ctrl_fillY, grayDark, 'FaceAlpha', 0.15, 'EdgeColor','none');

manip_fillX = [x_fit; flipud(x_fit)];
manip_fillY = [y_manip_ci(:,1); flipud(y_manip_ci(:,2))];
patch(manip_fillX, manip_fillY, orange, 'FaceAlpha', 0.15, 'EdgeColor','none');

% Fitted lines
plot(x_fit, y_ctrl_fit,  '-', 'Color', grayDark, 'LineWidth', 2);
plot(x_fit, y_manip_fit, '-', 'Color', orange,   'LineWidth', 2);

xlabel('Normalized mean firing rate');
ylabel('Normalized trajectory speed');
legend({'Baseline bins','eOPN bins','Baseline 95% CI','eOPN 95% CI', ...
        'Baseline fit','eOPN fit'}, 'Location','best');
title('Speed vs mean firing (session-averaged bins with linear fits + 95% CI)');
box off;
set(gca,'tickdir','out','fontsize',10);
axis square;

%%
% x_ctrl, y_ctrl  : baseline (Principal) binned points
% x_manip, y_manip: eOPN binned points

% Build one pooled table: FR, Condition, Speed
FR_all   = [x_ctrl(:);  x_manip(:)];          % column
Cond_all = [zeros(numel(x_ctrl),1); ...
            ones(numel(x_manip),1)];         % 0=baseline, 1=eOPN
Speed_all = [y_ctrl(:); y_manip(:)];

tbl_all = table(FR_all, Cond_all, Speed_all, ...
    'VariableNames', {'MeanFR','Condition','Speed'});

% Fit Speed ~ MeanFR + Condition on binned, session-averaged data
mdl_binned = fitlm(tbl_all, 'Speed ~ MeanFR + Condition');

disp(mdl_binned);


%%
% X: [N x 3], y: [N x 1], as you defined
% Design matrix: intercept, meanFR, condition
X = [ones(sum(valid),1), frValid, condValid];   % [N x 3]
y = spdValid;
% beta(2): effect of mean firing rate on speed
% beta(3): condition effect on speed after accounting for firing rate
[beta, betaCI, residuals, residualInt, stats] = regress(y, X);
% beta: [3 x 1]
% betaCI: [3 x 2] confidence intervals
% stats: [R2, F, p, errorVariance]
tbl = table(frValid, condValid, y, 'VariableNames', ...
    {'MeanFR','Condition','Speed'});

mdl = fitlm(tbl, 'Speed ~ MeanFR + Condition');
disp(mdl)
%% functions

function IntanBehaviour = grabTemp(IntanBehaviour)
for n = 1:IntanBehaviour.nCueHit
    IntanBehaviour.hitTemp(n,1) = IntanBehaviour.temperature(IntanBehaviour.cueHitTrace(n).LFPIndex(1));
end
IntanBehaviour.hitTemp = IntanBehaviour.hitTemp-IntanBehaviour.temperature(100);
for n = 1:IntanBehaviour.nCueMiss
    IntanBehaviour.missTemp(n,1) = IntanBehaviour.temperature(IntanBehaviour.cueMissTrace(n).LFPIndex(1));
end
IntanBehaviour.missTemp = IntanBehaviour.missTemp-IntanBehaviour.temperature(100);
for n = 1:length(IntanBehaviour.missTrace)
    IntanBehaviour.FATemp(n,1) = IntanBehaviour.temperature(IntanBehaviour.missTrace(n).LFPIndex(1));
end
IntanBehaviour.FATemp = IntanBehaviour.FATemp-IntanBehaviour.temperature(100);
end

% get mean matched spiking rates and trajectories
function matchedData = meanRateSpeed(Spikes,IntanBehaviour)
 % Number of trials and time bins
    nTrials = numel(Spikes.GPFA.hit.dat);
    [nNeurons, nTime] = size(Spikes.GPFA.hit.dat(1).spikes);
    [neuralDynamics,waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);

    % Mean firing rate per trial x time
    binSize = 20;  % 20 ms
    % Use the latent speed array to define nTime_binned
    [~, nTime_binned, ~] = size(neuralDynamics.hit.speed.speed);
    meanFR = nan(nTrials, nTime_binned);
    for tr = 1:nTrials
        S = Spikes.GPFA.hit.dat(tr).spikes;   % [nNeurons x nTime]
        [nNeurons, nTime1] = size(S);
        % Number of 20 ms bins
        nBins = floor(nTime1 / binSize);
        S = S(:, 1:(nBins * binSize));  % truncate to full bins

        % Reshape to [nNeurons x binSize x nBins] and sum across binSize
        S_reshaped = reshape(S, nNeurons, binSize, nBins);
        S_binned   = squeeze(sum(S_reshaped, 2));   % [nNeurons x nBins], spike counts per 20 ms
        meanFR(tr, :) = mean(S_binned, 1);          % average across neurons
    end

    % latentDim x time x trials
    spdLat = neuralDynamics.hit.speed.speed;

    % Combine across latent dimensions: sqrt(sum over dims of speed^2))
    trajSpeed = squeeze(sqrt(sum(spdLat.^2, 1)));   % [nTime x nTrials]
    trajSpeed = trajSpeed.';                        % [nTrials x nTime] to match meanFR

    % Define window in binned time indices
    tStart = 50;
    tEnd   = 120;

    % Crop to this window
    meanFR    = meanFR(:, tStart:tEnd);      % [nTrials x (tEnd-tStart+1)]
    trajSpeed = trajSpeed(:, tStart:tEnd);   % same size

    % Flatten to vectors: [nTrials*nTime x 1]
    frVec  = meanFR(:);
    spdVec = trajSpeed(:);

    % Find opto trials
    [IntanBehaviourBaseline,IntanBehaviourOpto, Waves, WavesOpto] = separateOptoTrials(IntanBehaviour,IntanBehaviour.parameters);
    baselineId = 1:length(IntanBehaviourBaseline.cueHitTrace);
    eOPNId = length(IntanBehaviourBaseline.cueHitTrace)+1:length(IntanBehaviour.cueHitTrace);
    assert(length(eOPNId)==length(IntanBehaviourOpto.cueHitTrace))

    % Trial-level condition labels
    % e.g., condLabel = zeros(nTrials,1); condLabel(thalTrials) = 1;
    condLabel = zeros(nTrials,1); 
    condLabel(baselineId) = 1;
    % Condition vector must match new time length
    nTime_win = size(meanFR, 2);
    condPerTrial = condLabel(:);
    condVec = repelem(condPerTrial, nTime_win, 1);

    % Define firing-rate bins (you can change # of bins)
    valid = ~isnan(frVec) & ~isnan(spdVec);
    frValid  = frVec(valid);
    spdValid = spdVec(valid);
    condValid = condVec(valid);

    matchedData.valid = valid;
    matchedData.frValid = frValid;
    matchedData.condValid = condValid;
    matchedData.spdValid = spdValid;
end