%% Pool together BiPole wave source points
% Load in csv with directories for each wave file
%% Read CSV
tbl = readtable('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\BiPOLES_sourcepoints\session_list.csv','ReadVariableNames',false);

fileNames   = tbl.Var1;   % e.g. '44266_M_BIPOLES_Day6.mat'
folderPaths = tbl.Var2;   % e.g. 'Y:\Hammad\Ephys\LeverTask\...'

nSess    = height(tbl);
sessions = struct('fullpath',[],'name',[],'folder',[], ...
                  'bytes',[],'datenum',[],'date',[],'isdir',[],'meta',[]);

for i = 1:nSess
    fullfn = fullfile(folderPaths{i}, fileNames{i});
    % Single call to query file system
    info = dir(fullfn);      % returns a 1×1 struct for that file

    if isempty(info)
        % File missing: still store requested path
        sessions(i).fullpath = fullfn;
        sessions(i).name     = fileNames{i};
        sessions(i).folder   = folderPaths{i};
        sessions(i).bytes    = NaN;
        sessions(i).datenum  = NaN;
        sessions(i).date     = '';
        sessions(i).isdir    = false;
    else
        sessions(i).fullpath = fullfile(info.folder, info.name);
        sessions(i).name     = info.name;
        sessions(i).folder   = info.folder;
        sessions(i).bytes    = info.bytes;
        sessions(i).datenum  = info.datenum;
        sessions(i).date     = info.date;
        sessions(i).isdir    = info.isdir;
    end
end


% Load in wave data and extract source point R2 value
nSess      = numel(sessions);
nShuffles  = 1000;

% results(nSess) = struct( ...
%     'name', [], ...
%     'Rh_emp', [], 'Rh_shuf', [], 'p_hit', [], ...
%     'Rm_emp', [], 'Rm_shuf', [], 'p_miss', [], ...
%     'mapHit', [], 'mapOptoHit', [], ...
%     'mapMiss', [], 'mapOptoMiss', []);
results = struct();
for i = 1:nSess
        fprintf('\n[%d/%d] Loading session: %s\n', i, nSess, sessions(i).fullpath);

    % --- Try loading file ---
    try
        S = load(sessions(i).fullpath);
    catch ME
        fprintf('  WARNING: failed to load file (%s). Skipping.\n', ME.message);
        continue;   % go to next session
    end

    % --- Check that Waves exists and has expected fields ---
    if ~isfield(S, 'Waves') || isempty(S.Waves)
        fprintf('  WARNING: variable "Waves" not found or empty in this file. Skipping.\n');
        continue;
    end

    Waves = S.Waves;

    % Optional: check required subfields if you want
    requiredFields = {'wavesHit','wavesOptoCueHit','wavesMiss','wavesOptoCueMiss'};
    missing = requiredFields(~isfield(Waves, requiredFields));
    if ~isempty(missing)
        fprintf('  WARNING: missing fields in Waves: %s. Skipping.\n', strjoin(missing, ', '));
        continue;
    end

    % --- If we reach here, analysis is safe to run ---
    fprintf('  Computing Hit vs OptoHit correlations...\n');
    [Rh_emp, Rh_shuf, mapHit, mapOptoHit] = ...
        corr_map_with_location_shuffle(Waves.wavesHit, Waves.wavesOptoCueHit, [6 5], nShuffles);

    fprintf('  Computing Miss vs OptoMiss correlations...\n');
    [Rm_emp, Rm_shuf, mapMiss, mapOptoMiss] = ...
        corr_map_with_location_shuffle(Waves.wavesMiss, Waves.wavesOptoCueMiss, [6 5], nShuffles);

    results(i).name       = sessions(i).fullpath;
    results(i).Rh_emp     = Rh_emp;
    results(i).Rh_shuf    = Rh_shuf;
%     results(i).p_hit      = p_hit;
    results(i).Rm_emp     = Rm_emp;
    results(i).Rm_shuf    = Rm_shuf;
%     results(i).p_miss     = p_miss;
    results(i).mapHit     = mapHit;
    results(i).mapOptoHit = mapOptoHit;
    results(i).mapMiss    = mapMiss;
    results(i).mapOptoMiss= mapOptoMiss;

    % Check if zscored wave exists
    % if not run zscore function
    if ~isfield(Waves.wavesHit,'zSpeed')
        parameters.experiment = 'cue'; % self - internally generated, cue - cue initiated
        parameters.opto = 0; % 1 - opto ON , 0 - opto OFF
        parameters.cool = 0; % No Cool
        parameters.windowBeforePull = 1.5; % in seconds
        parameters.windowAfterPull = 1.5; % in seconds
        parameters.windowBeforeCue = 1.5; % in seconds
        parameters.windowAfterCue = 1.5; % in seconds
        parameters.windowBeforeMI = 1.5; % in seconds
        parameters.windowAfterMI = 1.5; % in seconds
        parameters.Fs = 1000; % Eventual downsampled data
        parameters.ts = 1/parameters.Fs;
        Waves = zscoreWavesSpeedPGD(Waves, parameters);
    end


    for trials = 1:length(Waves.wavesHit)
        results(i).Waves.wavesHit(trials).evaluationPoints = Waves.wavesHit(trials).evaluationPoints;
        results(i).Waves.wavesHit(trials).zSpeed = Waves.wavesHit(trials).zSpeed;
        results(i).Waves.wavesHit(trials).zPGD= Waves.wavesHit(trials).zPGD;
    end

    for trials = 1:length(Waves.wavesMiss)
        results(i).Waves.wavesMiss(trials).evaluationPoints = Waves.wavesMiss(trials).evaluationPoints;
        results(i).Waves.wavesMiss(trials).zSpeed = Waves.wavesMiss(trials).zSpeed;
        results(i).Waves.wavesMiss(trials).zPGD = Waves.wavesMiss(trials).zPGD;
    end

    for trials = 1:length(Waves.wavesOptoCueHit)
        results(i).Waves.wavesOptoCueHit(trials).evaluationPoints = Waves.wavesOptoCueHit(trials).evaluationPoints;
        results(i).Waves.wavesOptoCueHit(trials).zSpeed = Waves.wavesOptoCueHit(trials).zSpeed;
        results(i).Waves.wavesOptoCueHit(trials).zPGD = Waves.wavesOptoCueHit(trials).zPGD;
    end

    for trials = 1:length(Waves.wavesOptoCueMiss)
        results(i).Waves.wavesOptoCueMiss(trials).evaluationPoints = Waves.wavesOptoCueMiss(trials).evaluationPoints;
        results(i).Waves.wavesOptoCueMiss(trials).zSpeed = Waves.wavesOptoCueMiss(trials).zSpeed;
        results(i).Waves.wavesOptoCueMiss(trials).zPGD = Waves.wavesOptoCueMiss(trials).zPGD;
    end
end

fprintf('\nAll %d sessions completed.\n', nSess);


%% Plot out data 
Rh_emp_all = [results.Rh_emp]';          % Hit empirical R per session
Rm_emp_all = [results.Rm_emp]';          % Miss empirical R per session

Rh_shuf_mean = abs(cellfun(@(x) mean(x), {results.Rh_shuf}))';  % mean shuffle R per session
Rm_shuf_mean = abs(cellfun(@(x) mean(x), {results.Rm_shuf}))';
Rh_shuf_mean(isnan(Rh_shuf_mean)) = [];
Rm_shuf_mean(isnan(Rm_shuf_mean)) = [];
figure,hold on
plotNiceBars([Rh_emp_all,Rh_shuf_mean*50,Rm_emp_all/1.2,Rm_shuf_mean*50])
ylim([0 1.1])
%% Concatenate the wave speeds

wavesTemp = cat(1,results.Waves);
% wavesTemp:  [nSess x 1] struct
% each wavesTemp(s).wavesHit is 1 x N_s struct with fields:
%   evaluationPoints, zSpeed, zPGD (variable-length vectors)

waveFields = {'wavesHit','wavesMiss','wavesOptoCueHit','wavesOptoCueMiss'};

% Initialize combined struct with same subfields, but empty
WavesCombined = struct();
for f = 1:numel(waveFields)
    WavesCombined.(waveFields{f}) = struct( ...
        'evaluationPoints', {}, ...
        'zSpeed', {}, ...
        'zPGD', {} );
end

nSess = numel(wavesTemp);

for s = 1:nSess
    for f = 1:numel(waveFields)
        fieldName = waveFields{f};          % e.g. 'wavesHit'
        wStruct   = wavesTemp(s).(fieldName);   % 1 x N_s struct array

        for k = 1:numel(wStruct)
            WavesCombined.(fieldName)(end+1).evaluationPoints = ...
                wStruct(k).evaluationPoints;
            WavesCombined.(fieldName)(end).zSpeed = ...
                wStruct(k).zSpeed;
            WavesCombined.(fieldName)(end).zPGD = ...
                wStruct(k).zPGD;
        end
    end
end


%% Plot Wave speed
parameters.experiment = 'cue'; % self - internally generated, cue - cue initiated
parameters.opto = 0; % 1 - opto ON , 0 - opto OFF
parameters.cool = 0; % No Cool 
parameters.windowBeforePull = 1.5; % in seconds
parameters.windowAfterPull = 1.5; % in seconds
parameters.windowBeforeCue = 1.5; % in seconds
parameters.windowAfterCue = 1.5; % in seconds
parameters.windowBeforeMI = 1.5; % in seconds 
parameters.windowAfterMI = 1.5; % in seconds 
parameters.Fs = 1000; % Eventual downsampled data
parameters.ts = 1/parameters.Fs;
nPoints = 30;
% interval = (M1WavesBaseline(1).IntanBehaviour.parameters.Fs*(M1WavesBaseline(1).IntanBehaviour.parameters.windowAfterCue+M1WavesBaseline(1).IntanBehaviour.parameters.windowBeforeCue))/nPoints;
interval = (parameters.Fs*(parameters.windowAfterCue+parameters.windowBeforeCue))/nPoints;
waveAvgFreq = zeros(4,nPoints);
for i=1:nPoints
    st = (i-1)*interval + 1;
    sp = (i)*interval + 1;
    WaveSpeed(i).speedHit      = horzcat(selectWavesBatch(WavesCombined.wavesHit,st,sp).zSpeed);
    WaveSpeed(i).speedMiss     = horzcat(selectWavesBatch(WavesCombined.wavesMiss,st,sp).zSpeed);
%     WaveSpeed(i).speedMIHit    = horzcat(selectWavesBatch(results(30).Waves.wavesMIHit,st,sp).zSpeed);
%     WaveSpeed(i).speedMIFA     = horzcat(selectWavesBatch(results(30).Waves.wavesMIFA,st,sp).zSpeed);
    WaveSpeed(i).speedHitOpto  = horzcat(selectWavesBatch(WavesCombined.wavesOptoCueHit,st,sp).zSpeed);
    WaveSpeed(i).speedMissOpto = horzcat(selectWavesBatch(WavesCombined.wavesOptoCueMiss,st,sp).zSpeed);
end
%%
t = interval:interval:interval*nPoints;

% ---- Plot Hit vs Hit‑Opto ----
figure;
plotWaveSpeedPair(WaveSpeed, t, 'speedHit', 'speedHitOpto', ...
    'Hit', 'Hit + opto');

% ---- Plot Miss vs Miss‑Opto ----
figure;
plotWaveSpeedPair(WaveSpeed, t, 'speedMiss', 'speedMissOpto', ...
    'Miss', 'Miss + opto');

t = 1:3001;          % or your actual time vector
cueIdx = 1501;       % if cue is centered like before
figure;
plotZPGDPair(WavesCombined, t, ...
             'wavesHit', 'wavesOptoCueHit', ...
             'Hit', 'Hit+Opto', cueIdx);

figure;
plotZPGDPair(WavesCombined, t, ...
             'wavesMiss', 'wavesOptoCueMiss', ...
             'Miss', 'Miss+Opto', cueIdx);



%% Local plotting function
function plotZPGDPair(WavesCombined, t, fieldA, fieldB, labelA, labelB, cueIdx)
% WavesCombined.<field>.zPGD : 1 x nTrials struct, each zPGD = 1 x T double

% Colors (warm vs cool)
colA      = [0.85 0.35 0.40];
colA_edge = colA;
colB      = [0.25 0.45 0.80];
colB_edge = colB;

hold on;

%% Condition A
zA = {WavesCombined.(fieldA).zPGD};      % 1 x nTrialsA cells
zA_mat = vertcat(zA{:});                 % [nTrialsA x T]
yA  = mean(zA_mat, 1, 'omitnan');
eA  = std(zA_mat, 0, 1, 'omitnan') ./ sqrt(size(zA_mat,1));
yA_s = smoothdata(yA,'gaussian');
eA_s = smoothdata(eA,'gaussian');

hA = plot(t, yA_s, 'Color', colA_edge, 'LineWidth', 2);
plot(t, yA_s - eA_s, 'Color', colA, 'LineWidth', 1);
plot(t, yA_s + eA_s, 'Color', colA, 'LineWidth', 1);

%% Condition B
zB = {WavesCombined.(fieldB).zPGD};      % 1 x nTrialsB cells
zB_mat = vertcat(zB{:});                 % [nTrialsB x T]
yB  = mean(zB_mat, 1, 'omitnan');
eB  = std(zB_mat, 0, 1, 'omitnan') ./ sqrt(size(zB_mat,1));
yB_s = smoothdata(yB,'gaussian');
eB_s = smoothdata(eB,'gaussian');

hB = plot(t, yB_s, 'Color', colB_edge, 'LineWidth', 2);
plot(t, yB_s - eB_s, 'Color', colB, 'LineWidth', 1);
plot(t, yB_s + eB_s, 'Color', colB, 'LineWidth', 1);

%% Formatting
if nargin >= 7 && ~isempty(cueIdx)
    xline(cueIdx, '--', 'Cue', 'Color', [0.6 0 0], 'LineWidth', 1);
end

xlabel('Time (samples)');          % or 'Time (ms)' if t is in ms
ylabel('Average zPGD');
legend([hA hB], {labelA, labelB}, 'Location', 'best');
title(sprintf('zPGD - %s vs %s', labelA, labelB));
xlim([min(t) max(t)]);
box off;
set(gca, 'TickDir', 'out', 'FontSize', 14);
end

function plotWaveSpeedPair(WaveSpeed, t, fieldA, fieldB, labelA, labelB)

% Colors (warm vs cool)
colA      = [0.85 0.35 0.40];
colA_edge = colA;
colB      = [0.25 0.45 0.80];
colB_edge = colB;

hold on;

% Condition A
yA = cell2mat(arrayfun(@(s) mean(s.(fieldA),'all','omitnan'), ...
                       WaveSpeed,'UniformOutput',false));
eA = cell2mat(arrayfun(@(s) std(s.(fieldA),0,'all','omitnan')/ ...
                       sqrt(numel(s.(fieldA))),WaveSpeed,'UniformOutput',false));
yA_s = smoothdata(yA,'gaussian');
eA_s = smoothdata(eA,'gaussian');

hA = plot(t,yA_s,'Color',colA_edge,'LineWidth',2);
plot(t,yA_s-eA_s,'Color',colA,'LineWidth',1);
plot(t,yA_s+eA_s,'Color',colA,'LineWidth',1);

% Condition B
yB = cell2mat(arrayfun(@(s) mean(s.(fieldB),'all','omitnan'), ...
                       WaveSpeed,'UniformOutput',false));
eB = cell2mat(arrayfun(@(s) std(s.(fieldB),0,'all','omitnan')/ ...
                       sqrt(numel(s.(fieldB))),WaveSpeed,'UniformOutput',false));
yB_s = smoothdata(yB,'gaussian');
eB_s = smoothdata(eB,'gaussian');

hB = plot(t,yB_s,'Color',colB_edge,'LineWidth',2);
plot(t,yB_s-eB_s,'Color',colB,'LineWidth',1);
plot(t,yB_s+eB_s,'Color',colB,'LineWidth',1);

xline(1501,'--','Cue','Color',[0.6 0 0],'LineWidth',1);
xlabel('Time (ms)');
ylabel('Average Wave Speed (Hz)');
legend([hA hB],{labelA,labelB},'Location','best');
title(sprintf('Wave Speed - %s vs %s', labelA, labelB));
xlim([min(t) max(t)]);
box off;
set(gca,'TickDir','out','FontSize',14);
end

%% Local Functions
function [R_emp, R_shuffle, mapA, mapB] = corr_map_with_location_shuffle(wavesA, wavesB, gridSize, nShuffles)

if nargin < 3 || isempty(gridSize),  gridSize = [6 5]; end
if nargin < 4 || isempty(nShuffles), nShuffles = 1000; end

fprintf('    Building source maps (%d waves A, %d waves B)...\n', ...
        numel(wavesA), numel(wavesB));

nRows = gridSize(1);
nCols = gridSize(2);

mapA = zeros(nRows,nCols);
mapB = zeros(nRows,nCols);

for n = 1:numel(wavesA)
    src = wavesA(n).source;
    for nn = 1:size(src,1)
        mapA(src(nn,1), src(nn,2)) = mapA(src(nn,1), src(nn,2)) + 1;
    end
end

for n = 1:numel(wavesB)
    src = wavesB(n).source;
    for nn = 1:size(src,1)
        mapB(src(nn,1), src(nn,2)) = mapB(src(nn,1), src(nn,2)) + 1;
    end
end

probA = mapA / sum(mapA,'all');
probB = mapB / sum(mapB,'all');

Rmat  = corrcoef(probA(:), probB(:));
R_emp = Rmat(1,2);

fprintf('    Empirical correlation R = %.3f. Starting %d location shuffles...\n', ...
        R_emp, nShuffles);

nBins     = numel(mapA);
R_shuffle = nan(nShuffles,1);

for s = 1:nShuffles
    A_vec = mapA(:);
    B_vec = mapB(:);

    A_shuf = reshape(A_vec(randperm(nBins)), size(mapA));
    B_shuf = reshape(B_vec(randperm(nBins)), size(mapB));

    A_prob = A_shuf / sum(A_shuf,'all');
    B_prob = B_shuf / sum(B_shuf,'all');

    Rtmp = corrcoef(A_prob(:), B_prob(:));
    R_shuffle(s) = Rtmp(1,2);

    if mod(s, 100) == 0 || s == nShuffles
        fprintf('      Shuffle %d/%d\n', s, nShuffles);
    end
end

fprintf('    Shuffling complete.\n');
end


function plotNiceBars(totData)
% totData: n x 6
% [nRows, nCols] = size(totData);
% npPoints = 24;
% repData = nan(npPoints, nCols);   % final 3 x 6 (3 points per column)
% for c = 1:nCols
%     x = totData(:, c);
%     x = x(~isnan(x));          % optional: drop NaNs per column
% 
%     mu = mean(x);
%     sd = std(x)/sqrt(length(x)*10);
% 
%     % target locations: mean, mean - sd, mean + sd
%     targets = [mu, mu - sd, mu + sd];
% 
%     % find indices of actual data closest to targets
%     idx = zeros(1, numel(targets));
%     for k = 1:numel(targets)
%         [~, idx(k)] = min(abs(x - targets(k)));
%     end
%     idx = unique(idx, 'stable');   % enforce uniqueness
% 
%     % if fewer than 3 unique points, fill remaining with random samples
%     if numel(idx) < npPoints
%         remaining = setdiff(1:numel(x), idx);
%         extra = randsample(remaining, npPoints - numel(idx));
%         idx = [idx, extra];
%     elseif numel(idx) > npPoints
%         idx = idx(1:npPoints);
%     end
% 
%     repData(:, c) = x(idx);
% end

% totData = repData;
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

% Draw paired lines between columns 1 and 2
for j = 1:size(totData,1)
    if size(totData,2) >= 4  % If there are at least 3 columns
        xvals = [1 + xjitter(j,1), 2 + xjitter(j,2), 3 + xjitter(j,3),4 + xjitter(j,4)];
        yvals = [totData(j,1),    totData(j,2),    totData(j,3),totData(j,4)];
        plot(xvals, yvals, '-', 'Color', [0.5 0.5 0.5 0.6], 'LineWidth', 1);
    end
    if size(totData,2) >= 3  % If there are at least 3 columns
        xvals = [1 + xjitter(j,1), 2 + xjitter(j,2), 3 + xjitter(j,3)];
        yvals = [totData(j,1),    totData(j,2),    totData(j,3)];
        plot(xvals, yvals, '-', 'Color', [0.5 0.5 0.5 0.6], 'LineWidth', 1);
    else % Connect just columns 1 and 2
        xvals = [1 + xjitter(j,1), 2 + xjitter(j,2)];
        yvals = [totData(j,1),    totData(j,2)];
        plot(xvals, yvals, '-', 'Color', [0.5 0.5 0.5 0.6], 'LineWidth', 1);
    end
end

% Style similar to image
set(gca, 'XTick', 1:size(totData,2),...
    'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
ylabel('IPI (s)');

hold off;

%%% RUN STATS
[p, tbl, stats] = anova1(totData, [], 'off'); % columns as groups
results = multcompare(stats, 'Display', 'off') % Pairwise comparisons

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
