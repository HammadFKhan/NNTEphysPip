% Compare similarities of source points across trial types
%% ---- Hit pair: wavesHit vs wavesOptoCueHit ----
[Rh_emp, Rh_shuf] = corr_map_with_location_shuffle( ...
    Waves.wavesHit, Waves.wavesOptoCueHit, [6 5], 1000);

%% ---- Miss pair: wavesMiss vs wavesOptoCueMiss ----
[Rm_emp, Rm_shuf] = corr_map_with_location_shuffle( ...
    Waves.wavesMiss, Waves.wavesOptoCueMiss, [6 5], 1000);

%% ---- p-values (one-sided, empirical > shuffle) ----
p_hit  = mean(Rh_shuf >= Rh_emp);
p_miss = mean(Rm_shuf >= Rm_emp);

fprintf('Hit:  R = %.3f, mean shuffle = %.3f, p = %.4f\n', ...
        Rh_emp, mean(Rh_shuf), p_hit);
fprintf('Miss: R = %.3f, mean shuffle = %.3f, p = %.4f\n', ...
        Rm_emp, mean(Rm_shuf), p_miss);

%% ---- Plot distributions if you want ----
figure;
subplot(1,2,1);
histogram(Rh_shuf); hold on;
xline(Rh_emp,'k','LineWidth',2);
title('Hit vs OptoHit');

subplot(1,2,2);
histogram(Rm_shuf); hold on;
xline(Rm_emp,'k','LineWidth',2);
title('Miss vs OptoMiss');

%% ---- Summary bar plot of R and shuffled mean ----
R_emp_all   = [Rh_emp, Rm_emp/2];
R_shuf_mean = [mean(Rh_shuf), mean(Rm_shuf)];

figure;
barData = [R_emp_all;];   % 2 (pairs) x 2 (emp vs shuffle)
plotNiceBars(barData)
set(gca,'XTickLabel',{'Hit','Miss'});
ylabel('Correlation');

legend({'Empirical R','Mean shuffled R'}, 'Location','best');
box off;

%% Local functions
%% ---- Helper function (put at bottom of file or in separate .m) ----
function [R_emp, R_shuffle] = corr_map_with_location_shuffle(wavesA, wavesB, gridSize, nShuffles)

if nargin < 3 || isempty(gridSize),  gridSize = [6 5]; end
if nargin < 4 || isempty(nShuffles), nShuffles = 1000; end

nRows = gridSize(1);
nCols = gridSize(2);

% Count maps
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

% Probabilities
probA = mapA / sum(mapA,'all');
probB = mapB / sum(mapB,'all');

% Empirical correlation
Rmat  = corrcoef(probA(:), probB(:));
R_emp = Rmat(1,2);

% Location shuffle
nBins     = numel(mapA);
R_shuffle = nan(nShuffles,1);

for s = 1:nShuffles
    % shuffle locations within each map
    A_vec = mapA(:);
    B_vec = mapB(:);

    A_shuf = reshape(A_vec(randperm(nBins)), size(mapA));
    B_shuf = reshape(B_vec(randperm(nBins)), size(mapB));

    A_prob = A_shuf / sum(A_shuf,'all');
    B_prob = B_shuf / sum(B_shuf,'all');

    Rtmp = corrcoef(A_prob(:), B_prob(:));
    R_shuffle(s) = Rtmp(1,2);
end
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
totData  = totData(:,2:end);
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