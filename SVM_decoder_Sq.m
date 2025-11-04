% Build an SVM decoder for logistic decoding of task outcome
% Sort hit and effort trials
latentDynamics = cat(3,M1neuralDynamics.hiteffort.X,M1neuralDynamics.effort.X);
trialLabels = [ones(size(M1neuralDynamics.hiteffort.X,3),1);zeros(size(M1neuralDynamics.effort.X,3),1)];
[F1effort,F1effort_shuf] = getSVMdecode(latentDynamics,trialLabels);

latentDynamics = cat(3,M1neuralDynamics.hit.X,M1neuralDynamics.miss.X);
trialLabels = [ones(size(M1neuralDynamics.hit.X,3),1);zeros(size(M1neuralDynamics.miss.X,3),1)];
[F1miss,F1miss_shuf] = getSVMdecode(latentDynamics,trialLabels);
%%
figure; hold on
nShuffles = length(F1effort_shuf);
bar(1, F1effort, 0.5, 'FaceColor', [0.7 0.3 0.3]);
bar(2, mean(F1effort_shuf), 0.5, 'FaceColor', [0.7 0.7 0.7]);
% Overlay scatter points for shuffle scores
scatter(ones(nShuffles,1)*2, F1effort_shuf, 40, 'k', 'filled', 'jitter','on','jitterAmount',0.15);
% Add error bars
errorbar([1 2], [F1effort mean(F1effort_shuf)], ...
    [0 std(F1effort_shuf)/sqrt(nShuffles)], 'k', 'LineStyle', 'none');

xticks([1 2]);
xticklabels({'',''});
ylabel('Decoder performance (F1 Score)');
title('Hit vs Effort Decoder Performance');
% Optionally, run t-test and annotate
[~,p] = ttest(repmat(F1effort,nShuffles,1), F1effort_shuf);
p_str = sprintf('p = %.3g', p); % Formats p to 3 significant digits

if p < 0.05
    text(1.5, max([F1effort; F1effort_shuf]) + 0.02, ...
        ['* (' p_str ')'], 'HorizontalAlignment', 'center', 'FontSize', 12);
else
    text(1.5, max([F1effort; F1effort_shuf]) + 0.02, ...
        ['n.s. (' p_str ')'], 'HorizontalAlignment', 'center', 'FontSize', 12);
end
hold off
set(gca,'TickDir', 'out', 'Box', 'off', 'FontSize', 12),box off, axis square
%% Get Hit vs miss
figure; hold on

bar(1, F1miss, 0.5, 'FaceColor', [0.7 0.3 0.3]);
bar(2, mean(F1miss_shuf), 0.5, 'FaceColor', [0.7 0.7 0.7]);
% Overlay scatter points for shuffle scores
scatter(ones(nShuffles,1)*2, F1miss_shuf, 40, 'k', 'filled', 'jitter','on','jitterAmount',0.15);
% Add error bars
errorbar([1 2], [F1miss mean(F1miss_shuf)], ...
    [0 std(F1miss_shuf)/sqrt(nShuffles)], 'k', 'LineStyle', 'none');

xticks([1 2]);
xticklabels({'',''});
ylabel('Decoder performance (F1 Score)');
title('Hit vs Miss Decoder Performance');
% Optionally, run t-test and annotate
[~,p] = ttest(repmat(F1miss,nShuffles,1), F1miss_shuf);
p_str = sprintf('p = %.3g', p); % Formats p to 3 significant digits

if p < 0.05
    text(1.5, max([F1miss; F1miss_shuf]) + 0.02, ...
        ['* (' p_str ')'], 'HorizontalAlignment', 'center', 'FontSize', 12);
else
    text(1.5, max([F1miss; F1miss_shuf]) + 0.02, ...
        ['n.s. (' p_str ')'], 'HorizontalAlignment', 'center', 'FontSize', 12);
end
hold off
set(gca,'TickDir', 'out', 'Box', 'off', 'FontSize', 12),box off, axis square
%%
function [decodeScore_true,decodeScores_shuffle] = getSVMdecode(latentDynamics,trialLabels)
% Assume X is trials x features matrix, Y is trials x 1 vector of labels

[numDims, numTimepoints, numTrials] = size(latentDynamics);
X = reshape(latentDynamics, [numDims * numTimepoints, numTrials])'; % trials x features

Y = trialLabels; % trial outcome labels

% Calculate class frequencies
numEffort = sum(Y == 0);
numHit = sum(Y == 1);

% Inverse frequency weights
weightEffort = numHit / (numEffort + numHit);
weightHit = numEffort / (numEffort + numHit);

% Cost matrix for SVM (rows: true class, columns: predicted class)
costMatrix = [0, weightEffort; weightHit, 0];
% Step 6: Train a linear SVM classifier on the flattened trial features
SVMModel = fitcsvm(X, Y, 'KernelFunction', 'linear','Cost', costMatrix);

% Step 7: Perform k-fold cross-validation to estimate out-of-sample accuracy
CVSVMModel = crossval(SVMModel, 'Leaveout', 'on');
classLoss = kfoldLoss(CVSVMModel);

fprintf('Cross-validated classification loss (error rate): %.3f\n', classLoss);
fprintf('Cross-validated accuracy: %.3f\n', 1 - classLoss);

labelsPred = kfoldPredict(CVSVMModel);
labelsPred_train = predict(SVMModel, X); % In-sample prediction (may be inflated accuracy)

acc_heldout = mean(labelsPred == Y);

% Accuracy or F1 score for in-sample (trained) data
acc_train = mean(labelsPred_train == Y);

fprintf('Held-out (CV) accuracy: %.3f\n', acc_heldout);
fprintf('Trained (in-sample) accuracy: %.3f\n', acc_train);

% Calculate F1-score
% y_true: ground truth labels (e.g., 0/1)
% y_pred: model predictions (e.g., 0/1)
[decodeScore_true] = f1score(labelsPred, Y);
%% Shuffled response
nShuffles = 10;
decodeScores_shuffle = zeros(nShuffles,1);
[numDims, numTimepoints, numTrials] = size(latentDynamics);
X = reshape(latentDynamics, [numDims * numTimepoints, numTrials])'; % trials x features
for s = 1:nShuffles
    Y_shuf = Y(randperm(length(Y))); % Shuffle trial labels
    SVMModel_shuf = fitcsvm(X, Y_shuf, 'KernelFunction', 'linear','Cost', costMatrix);
    CVSVMModel_shuf = crossval(SVMModel_shuf);
    labelsPred_shuf = kfoldPredict(CVSVMModel_shuf);
    % F1 or accuracy
    decodeScores_shuffle(s) = f1score(labelsPred_shuf, Y_shuf);
end
end

function [f1, precision, recall] = f1score(y_true, y_pred)
    % Returns F1, precision, and recall for binary labels
    TP = sum((y_true == 1) & (y_pred == 1));
    FP = sum((y_true == 0) & (y_pred == 1));
    FN = sum((y_true == 1) & (y_pred == 0));

    if TP+FP == 0
        precision = 0;
    else
        precision = TP / (TP + FP);
    end
    if TP+FN == 0
        recall = 0;
    else
        recall = TP / (TP + FN);
    end
    if precision+recall == 0
        f1 = 0;
    else
        f1 = 2 * (precision * recall) / (precision + recall);
    end
end

function plotNiceBars(totData)
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
set(gca, 'XTick', 1:size(totData,2), 'XTickLabel', {'Second Pull', 'Third Pull', 'Polymer', 'Late'}, ...
    'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
ylabel('IPI (s)');
ylim([0 2]);

hold off;

%%% RUN STATS
[p, tbl, stats] = anova1(totData, [], 'off'); % columns as groups
results = multcompare(stats, 'Display', 'off'); % Pairwise comparisons

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
%% Now do the same thing for miss
