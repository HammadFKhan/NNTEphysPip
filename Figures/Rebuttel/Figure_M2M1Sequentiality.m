%%
% Calculate M2 and M1 diff
load('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\M1M2SequentialityCool\M1Sequentiality.mat')
load('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\M1M2SequentialityCool\M2Sequentiality.mat')
%%%
dat = [];
for n = 1:length(M2baseline)
    if isfield(M2baseline(n).baselineSqEntropy.CueHit,'SqI') && isfield(M2cooling(n).cooledSqEntropy.CueHit,'SqI')
        dat = vertcat(dat,[M2baseline(n).baselineSqEntropy.CueHit.SqI(1)',M2cooling(n).cooledSqEntropy.CueHit.SqI(3)']);
    end
end
M2SqI = diff(dat,1,2);
dat = [];
for n = 1:length(M1baseline)
    if isfield(M1baseline(n).baselineSqEntropy.CueHit,'SqI') && isfield(M1cooling(n).cooledSqEntropy.CueHit,'SqI')
        dat = vertcat(dat,[M1baseline(n).baselineSqEntropy.CueHit.SqI(1)',M1cooling(n).cooledSqEntropy.CueHit.SqI(3)']);
    end
end
M1SqI = diff(dat,1,2);
SqI_tot = nan(max([length(M1SqI), length(M2SqI)]),2);
SqI_tot(1:length(M1SqI),1) = M1SqI;
SqI_tot(1:length(M2SqI),2) = M2SqI;
plotDat(SqI_tot)
%% LOCAL FUNCTION
function plotDat(dat)
baseline_data = dat(:,1);
eOPN_data = dat(:,2);
figure;hold on % Create a new figure window

% Define colors for the points
baseline_color = [0.6 0.6 0.6]; % Gray
cooling_color = [0.1 0.35 0.85];   % Orange

% Plot lines connecting paired points first (light gray, thin)
% hold on; % Keep the plot active for multiple elements
% for i = 1:numel(baseline_data)
%     plot([1, 2], [baseline_data(i), eOPN_data(i)], 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5); % Light gray line
% end

% Plot individual data points using scatter
% Baseline points (x=1)
scatter(ones(size(baseline_data)), baseline_data, 100, 'filled', ...
        'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', baseline_color);

% eOPN points (x=2)
scatter(2 * ones(size(eOPN_data)), eOPN_data, 100, 'filled', ...
        'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', cooling_color);

% Add a horizontal dashed line at y = 0 for reference
%plot(xlim, [0 0], 'k--', 'LineWidth', 1.5);

% --- Beautify the plot ---
ax = gca; % Get current axes handle

ax.TickDir = 'out'; % Ticks point outwards
ax.FontSize = 14;   % Font size for tick labels
ax.Box = 'off';     % Turn off the box around the plot

% Set x-axis limits and labels for two groups
xlim([0.5 2.5]); % Adjust limits to center the two columns
xticks([1 2]); % Set tick marks at x=1 and x=2
xticklabels({'M1', 'M2'}); % Set x-axis labels
ax.XAxis.Color = [0.3 0.3 0.3]; % Darker color for x-axis labels
ax.YAxis.Color = [0.3 0.3 0.3]; % Darker color for y-axis labels
xlabel(''); % No overall x-axis label needed

% Y-axis label (match the image more closely)
ylabel('SI Index', 'FontSize', 16); 

% ylim([0.6 0.9]); % Adjust limits to ensure 0 is visible

% --- Add statistical annotation (line and p-value) ---
% Get current y-axis limits to place the p-value
yLimits = ylim(ax);
xLimits = xlim(ax);

% Position for the p-value line and text (adjust these values manually for best fit)
y_line = yLimits(2) * 0.95; % Near the top
x_left = 1; % Corresponds to Baseline x-position
x_right = 2; % Corresponds to eOPN x-position

% Draw the line connecting the two groups for annotation
line([x_left, x_right], [y_line, y_line], 'Color', 'k', 'LineWidth', 0.5);

% Draw the small vertical bars at the ends of the horizontal line
line([x_left, x_left], [y_line - (range(yLimits)*0.02), y_line], 'Color', 'k', 'LineWidth', 0.5);
line([x_right, x_right], [y_line - (range(yLimits)*0.02), y_line], 'Color', 'k', 'LineWidth', 0.5);

% Add the p-value text (using the paired t-test p-value as in your second image)
text_x_pos = (x_left + x_right) / 2;
text_y_pos = y_line + (range(yLimits)*0.03); % Slightly above the line
[h_ttest_paired, p_ttest_paired, ci_ttest_paired, stats_ttest_paired] = ttest2(eOPN_data, baseline_data);
text(text_x_pos, text_y_pos, sprintf('p = %.4f', p_ttest_paired), ...
     'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'bottom', ...
     'FontSize', 14, 'FontWeight', 'bold'); % Similar font size/weight as image

% Optional: Add a title if you want, but the image you provided doesn't have one
% title('TW Speed Modulation', 'FontSize', 16);

hold off; % Release the plot
end

