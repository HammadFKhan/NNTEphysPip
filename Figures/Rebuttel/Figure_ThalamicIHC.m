%% IHC analysis of thalamic nuclei

%% Load table
fileName = "Y:\Hammad\Ephys\LeverTask\LeverTaskRebuttal\IHCThalamus\PooledData.xlsx";
T = readtable(fileName);

% Convert selected columns to matrix
% Example: use columns 2–4 of the table (edit these to your nuclei columns)
colsUse = [1:4,6];
rawData = T{:, colsUse};          % numeric matrix [nAnimals x nConditions]
rawData(:,[1,5]) = rawData(:,[1,5])+2314;
% Row-normalize (each row sums to 1)
rowSums = nansum(rawData, 2);
normData = rawData ./ rowSums;    % same size as rawData

%% Plot normalized data
figure; 
plotNiceBars_noLines(normData);   % custom function below
set(gca, 'XTickLabel', T.Properties.VariableNames(colsUse)); % label nuclei
ylabel('Normalized fluorescence');
axis square
ylim([0 0.5])
%% LOCAL FUNCTION
function plotNiceBars_noLines(totData)

if istable(totData)
    totData = totData{:,:};
end

means = nanmean(totData);                         % Bar heights
sems  = nanstd(totData) ./ sqrt(size(totData,1)); % SEM

hold on;

% Bars
b = bar(means, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k');

% Error bars
errorbar(1:size(totData,2), means, sems, 'k', ...
    'LineStyle', 'none', 'LineWidth', 1);

% Jittered points
xjitter = randn(size(totData))*0.05;
for i = 1:size(totData,2)
    scatter(i + xjitter(:,i), totData(:,i), 25, 'o', ...
        'MarkerEdgeColor', [0.25 0.25 0.25], ...
        'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4);
end

set(gca, 'XTick', 1:size(totData,2), ...
    'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

ylabel('Value');

%%% Stats: one-way ANOVA across columns
[p, tbl, stats] = anova1(totData, [], 'off');   % columns as groups
results = multcompare(stats, 'Display', 'off'); % pairwise

disp(['ANOVA p-value: ', num2str(p)]);
alpha = 0.05;
sigPairs = results(results(:,6) < alpha, :);

% Significance bars
baseY      = max(means + sems) * 1.05;
offsetStep = max(means + sems) * 0.05;

if isempty(sigPairs)
    % no significant pairwise effects: report omnibus F, p
    Fstat   = cell2mat(tbl(2,5));
    p_anova = p;
    xPos = size(totData,2)/2;
    yPos = max(means + sems) * 1.3;
    text(xPos, yPos, sprintf('ANOVA F=%.2f, p=%.3f', Fstat, p_anova), ...
        'HorizontalAlignment', 'center', 'FontSize', 10);
else
    for i = 1:size(sigPairs,1)
        x1 = sigPairs(i,1);
        x2 = sigPairs(i,2);
        y  = baseY + (i-1)*offsetStep;

        plot([x1 x1 x2 x2], [y y+offsetStep y+offsetStep y], ...
             'k-', 'LineWidth', 1);

        text(mean([x1 x2]), y + offsetStep*0.1, '*', ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 16, 'FontWeight', 'bold');
    end
    xPos = size(totData,2)/2;
    yPos = max(means + sems) * 1.3;
    Fstat   = cell2mat(tbl(2,5));
    p_anova = p;
    text(xPos, yPos, sprintf('ANOVA F=%.2f, p=%.3f', Fstat, p_anova), ...
        'HorizontalAlignment', 'center', 'FontSize', 10);
end

hold off;

end