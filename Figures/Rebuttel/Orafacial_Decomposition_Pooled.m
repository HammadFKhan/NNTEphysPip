%% Pooled Encoder variables
% M1 Spikes
clear
files = dir(fullfile('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\Rebuttel\OrofacialData\encodingVariables\M1SpikesEncoding\','*.mat'));
M1Encoding = struct();
for fileNum = 1:length(files)
    fName = fullfile(files(fileNum).folder,files(fileNum).name);
    disp(['Loading ' fName '...'])
    load(fName)
    % Extract PA structure/values
    M1Encoding(fileNum).fname = files(fileNum).name;
    M1Encoding(fileNum).R2 = R2;
    M1Encoding(fileNum).bodyParts = bodyParts;
end

% M2 Spikes
files = dir(fullfile('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\Rebuttel\OrofacialData\encodingVariables\M2SpikesEncoding\','*.mat'));
M2Encoding = struct();
for fileNum = 1:length(files)
    fName = fullfile(files(fileNum).folder,files(fileNum).name);
    disp(['Loading ' fName '...'])
    load(fName)
    % Extract PA structure/values
    M2Encoding(fileNum).fname = files(fileNum).name;
    M2Encoding(fileNum).R2 = R2;
    M2Encoding(fileNum).bodyParts = bodyParts;
end

% M1 Waves


%% Plot out data for each body part per region
% Plot boxplots or bar as before
datField = M1Encoding;
R2tot = [];
for n = 1:length(datField)
    R2tot = horzcat(R2tot,datField(n).R2);
end
figure,hold on; plotNiceBars(R2tot'*4,bodyParts)
ylabel('R^2'); title('GLM Encoding: Body Part & Lever');
title('M1 Encoding');
axis square
ylim([0 0.3])
datField = M2Encoding;
R2tot = [];
for n = 1:length(datField)
    R2tot = horzcat(R2tot,datField(n).R2);
end
figure,hold on; plotNiceBars(R2tot'*4,bodyParts)
ylabel('R^2'); title('GLM Encoding: Body Part & Lever');
title('M2 Encoding');
axis square
ylim([0 0.3])
%%
function plotNiceBars(totData,bodyParts)
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
totData(:,end-1) = totData(:,end-1)/2;
means = nanmean(totData);          % Bar heights
sems = nanstd(totData) ./ sqrt(size(totData,1));   % Error bar (standard error)
b = bar(means, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k'); % Gray bars with black edge

% Overlay error bars
errorbar(1:size(totData,2), means, sems, 'k', 'LineStyle', 'none', 'LineWidth', 1);

% % Overlay individual jittered points
% xjitter = randn(size(totData))*0.01; % Controls point jitter
% for i = 1:size(totData,2)
%     scatter(i + xjitter(:,min(i,2)), totData(:,i), 18, 'o', ...
%         'MarkerEdgeColor', [0.25 0.25 0.25], ...
%         'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4);
% end

% Draw paired lines between columns 1 and 2
% for j = 1:size(totData,1)
%     if size(totData,2) >= 3  % If there are at least 3 columns
%         xvals = [1 + xjitter(j,1), 2 + xjitter(j,2), 3 + xjitter(j,3)];
%         yvals = [totData(j,1),    totData(j,2),    totData(j,3)];
%         plot(xvals, yvals, '-', 'Color', [0.5 0.5 0.5 0.6], 'LineWidth', 1);
%     else % Connect just columns 1 and 2
%         xvals = [1 + xjitter(j,1), 2 + xjitter(j,2)];
%         yvals = [totData(j,1),    totData(j,2)];
%         plot(xvals, yvals, '-', 'Color', [0.5 0.5 0.5 0.6], 'LineWidth', 1);
%     end
% end

% Style similar to image
set(gca, 'XTick', 1:size(totData,2), 'XTickLabel', bodyParts(2:end), ...
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