%% Intrinsic manifold analysis
tempCutoff = -12; % Cuttoff of temperature cooling
hitbaseline = [];hitcooled = [];
missbaseline = [];misscooled = [];
FAbaseline = [];FAcooled = [];
dimension = 1;
for n = 1:length(M1DynamicsCooled)
    temperatureId = (M1DynamicsCooled(n).IntanBehaviour.hitTemp<tempCutoff);
    hitbaseline{n} = M1DynamicsCooled(n).neuralDynamics.hit.s(dimension,~temperatureId);
    hitcooled{n} = M1DynamicsCooled(n).neuralDynamics.hit.s(dimension,temperatureId);
end
hitbaseline = cellfun(@mean,hitbaseline);
hitcooled = cellfun(@mean,hitcooled);
% hitcooled(end) = hitbaseline(end)-0.1;
% hitcooled(5) = hitbaseline(5)+0.6;

%%
totData = [hitbaseline',hitcooled']+0.3;
figure,hold on
plotNiceBars(totData)
ylim([0 0.8])
%%
figure,customBarplot([hitbaseline',hitcooled']);
box off,set(gca,'tickdir','out','fontsize',14),axis square
ylabel('Intrinsic Trajectory Manifold')
ylim([0 10])
axis square
[h,p] = ttest2(hitbaseline,hitcooled);

function customBarplot(data,varargin)
if ~isempty(varargin) && (strcmp(varargin{1},'Scatter') || strcmp(varargin{1},'scatter'))
    if strcmp(varargin{2},'on')
        scatterOn = 1;
    else
        scatterOn = 0;
    end
else
    scatterOn = 1;
end

labels = []; buff = [];
if scatterOn
    hold on
    % Plot scatter points with jitter
    for i = 1:size(data,2)
        t = data(:,i);
        scatter(i*ones(length(t(t~=0)),1), t(t~=0), 'filled',...
            'jitter','on','jitterAmount',0.0), hold on
    end
    
    % Add connecting lines for each row
    for row = 1:size(data,1)
        x = [];
        y = [];
        for col = 1:size(data,2)
            if data(row,col) ~= 0
                x(end+1) = col;
                y(end+1) = data(row,col);
            end
        end
        if length(x) > 1
            plot(x, y, '-k', 'LineWidth', 0.5) % Connect points with black lines
        end
    end
end

% Plot bars and errorbars
for i = 1:size(data,2)
    t = data(:,i);
    buff = t(t~=0);
    labels = i*ones(length(buff),1);
    bar(i, nanmean(buff)), hold on
    err = nanstd(buff)/sqrt(length(buff));
    errorbar(i, nanmean(buff), err), hold on
end

h = findobj('LineStyle','--'); 
set(h, 'LineStyle','-');
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

nTrials = size(totData,1);
% Draw paired lines between columns 1 and 2
for j = 1:size(totData,1)
    if size(totData,2) >= 3  % If there are at least 3 columns
        xvals = [1 + xjitter(j,1), 2 + xjitter(j,2), 3 + xjitter(j,3)];
        yvals = [totData(j,1),    totData(j,2),    totData(j,3)];
        plot(xvals, yvals, '-', 'Color',[0.3 0.3 0.3],'LineWidth', 1);
    else % Connect just columns 1 and 2
        xvals = [1 + xjitter(j,1), 2 + xjitter(j,2)];
        yvals = [totData(j,1),    totData(j,2)];
        plot(xvals, yvals, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1);
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

