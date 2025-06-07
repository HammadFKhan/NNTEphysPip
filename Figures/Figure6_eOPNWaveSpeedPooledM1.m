% 1. Define your data
% Assuming 'modulation_data' represents the 'eOPN' state or modulated state
% and 'baseline_data' represents the 'Baseline' state, for the same 5 samples.
% If 'modulation_data' is ALREADY the difference, then we treat baseline_data
% as just raw values to plot for comparison. Let's assume they are two separate
% measurements from the same subjects/samples.
load('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\eOPN\M1WaveSpeedPooled.mat')

% 2. Perform Paired Statistical Tests
% We are now testing the difference between the 'eOPN_data' and 'baseline_data'
% against a null hypothesis of zero difference.

% --- Wilcoxon Signed-Rank Test (Paired) ---
% Compares two related (paired) samples. It checks if the median of the differences
% between paired observations is significantly different from zero.
[p_wilcoxon_paired, h_wilcoxon_paired, stats_wilcoxon_paired] = signrank(eOPN_data, baseline_data);

% --- Paired Sample t-test ---
% Compares the means of two related (paired) samples. It assumes the differences
% between pairs are normally distributed.
[h_ttest_paired, p_ttest_paired, ci_ttest_paired, stats_ttest_paired] = ttest(eOPN_data, baseline_data);

% 3. Plot the data points as two vertical columns

figure; % Create a new figure

% Plot Baseline data (left column, e.g., x=1)
scatter(ones(size(baseline_data)), baseline_data, 100, 'filled', ...
        'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', [0.6 0.6 0.6]); % Gray for Baseline
hold on; % Keep the plot active to add more elements

% Plot eOPN data (right column, e.g., x=2)
scatter(2 * ones(size(eOPN_data)), eOPN_data, 100, 'filled', ...
        'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', [0.85 0.35 0.1]); % Orange for eOPN

% Optionally, plot lines connecting paired points
% This is especially good for paired data to visualize individual changes
for i = 1:numel(baseline_data)
    plot([1, 2], [baseline_data(i), eOPN_data(i)], 'k-', 'LineWidth', 0.5, 'Color', [0.7 0.7 0.7]); % Light gray line
end

% Re-plot the scatter points on top of the lines to ensure visibility
scatter(ones(size(baseline_data)), baseline_data, 100, 'filled', ...
        'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', [0.6 0.6 0.6]); % Gray for Baseline
scatter(2 * ones(size(eOPN_data)), eOPN_data, 100, 'filled', ...
        'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', [0.85 0.35 0.1]); % Orange for eOPN


% Add a horizontal line at y = 0 for reference
plot(xlim, [0 0], 'k--', 'LineWidth', 1.5); % Black dashed line

% --- Beautify the plot ---
ax = gca; % Get current axes handle

ax.TickDir = 'out'; % Ticks point outwards
ax.FontSize = 14;   % Font size for tick labels
ax.Box = 'off';     % Turn off the box around the plot

% Set x-axis limits and labels for two groups
xlim([0.5 2.5]); % Adjust limits to center the two columns
xticks([1 2]); % Set tick marks at x=1 and x=2
xticklabels({'Baseline', 'eOPN'}); % Set x-axis labels
xlabel(''); % No overall x-axis label needed

% Y-axis label
ylabel('Speed modulation', 'FontSize', 16); % Change 'Units' to whatever your data represents

% Set y-axis limits to clearly show all points and the zero line
all_data = [baseline_data, eOPN_data];
min_val = min(all_data) - 0.1;
max_val = max(all_data) + 0.1;
if min_val > 0
    min_val = -0.1; % Ensure zero is visible even if all points are positive
end
if max_val < 0
    max_val = 0.1; % Ensure zero is visible even if all points are negative
end
ylim([min(0, min_val) max(0, max_val)]);

% Add a title
title('Comparison of Baseline vs. eOPN Values', 'FontSize', 16);

% Add the p-value annotation for the paired test
% Determine a good position for the p-value text above the columns
y_line_pos = max(ax.YLim) * 0.95; % Near the top of the plot
x_line_start = 1;
x_line_end = 2;

% Draw the line connecting the two groups for annotation
line([x_line_start, x_line_end], [y_line_pos, y_line_pos], 'Color', 'k', 'LineWidth', 0.5);
line([x_line_start, x_line_start], [y_line_pos - (range(ax.YLim)*0.02), y_line_pos], 'Color', 'k', 'LineWidth', 0.5);
line([x_line_end, x_line_end], [y_line_pos - (range(ax.YLim)*0.02), y_line_pos], 'Color', 'k', 'LineWidth', 0.5);


text_x_pos = (x_line_start + x_line_end) / 2;
text_y_pos = y_line_pos + (range(ax.YLim)*0.03);


% Optionally, add the paired t-test p-value as well
text(text_x_pos, text_y_pos - (range(ax.YLim)*0.05), ... % Slightly lower position
     sprintf('Paired t-test p = %.4f', p_ttest_paired), ...
     'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'bottom', ...
     'FontSize', 10, 'Color', [0.4 0.4 0.4]);

hold off; % Release the plot

% Display results in the command window
fprintf('--- Statistical Test Results (Null Hypothesis: No Difference Between Paired Samples) ---\n');
fprintf('Baseline Data: [%s]\n', num2str(baseline_data));
fprintf('eOPN Data: [%s]\n', num2str(eOPN_data));

fprintf('\nPaired Sample t-test:\n');
fprintf('  p-value: %.4f\n', p_ttest_paired);
fprintf('  Significance (h=1 if significant at 5%% level): %d\n', h_ttest_paired);
fprintf('  t-statistic (from stats struct): %.2f\n', stats_ttest_paired.tstat);
fprintf('  Degrees of freedom: %d\n', stats_ttest_paired.df);
axis square
ylim([-0.8 0.2 ])