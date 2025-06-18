%%
load("\\10.165.57.13\Sutter_backup\Hammad\Ephys\LeverTask\Data_for_Figures\WavesCooling\ReactionTime.mat")
% Assuming your data structure is named 'data' and has fields for each condition
% Adjust field names based on your actual structure
baseline_RT = [];
cooled_RT = [];
recovery_RT = [];
for n = 1:6
    % Extract reaction time data for each condition
    baseline_RT = horzcat(baseline_RT,RT{n}{1});  % or whatever your baseline field is named
    cooled_RT = horzcat(cooled_RT,RT{n}{2});       % or whatever your cooled field is named
    recovery_RT = horzcat(recovery_RT,RT{n}{3});   % or whatever your recovery field is named
end
% Calculate baseline mean
baseline_mean = mean(baseline_RT, 'omitnan');

% Normalize data as percent change from baseline
baseline_norm = ((baseline_RT - baseline_mean) / baseline_mean) * 100;
cooled_norm = ((cooled_RT - baseline_mean) / baseline_mean) * 100;
recovery_norm = ((recovery_RT - baseline_mean) / baseline_mean) * 100;

% Alternative: Z-score normalization
% baseline_std = std(baseline_RT, 'omitnan');
% baseline_norm = (baseline_RT - baseline_mean) / baseline_std;
% cooled_norm = (cooled_RT - baseline_mean) / baseline_std;
% recovery_norm = (recovery_RT - baseline_mean) / baseline_std;

% Prepare data for plotting
all_data = [baseline_norm'; cooled_norm'; recovery_norm'];
group_labels = [repmat({'Baseline'}, length(baseline_norm), 1); ...
                repmat({'Cooled'}, length(cooled_norm), 1); ...
                repmat({'Recovery'}, length(recovery_norm), 1)];

% Prepare data for analysis
conditions = {'Baseline', 'Cooled', 'Recovery'};
data_cell = {baseline_norm, cooled_norm, recovery_norm};
n_conditions = length(data_cell);

%%% 3. Prepare data for omnibus tests
% Create vectors for ANOVA/Kruskal-Wallis
all_data = [];
group_labels = [];
for i = 1:n_conditions
    all_data = [all_data; data_cell{i}(:)];
    group_labels = [group_labels; i*ones(length(data_cell{i}), 1)];
end
fprintf('=== OMNIBUS TESTS ===\n');

% One-way ANOVA (parametric)
[p_anova, tbl_anova, stats_anova] = anova1(all_data, group_labels, 'off');
fprintf('One-way ANOVA: F(%d,%d) = %.3f, p = %.4f\n', ...
        tbl_anova{2,3}, tbl_anova{4,3}, tbl_anova{2,5}, p_anova);

% Kruskal-Wallis test (non-parametric)
[p_kw, tbl_kw, stats_kw] = kruskalwallis(all_data, group_labels, 'off');
fprintf('Kruskal-Wallis: χ²(%d) = %.3f, p = %.4f\n', ...
        tbl_kw{2,3}, tbl_kw{2,5}, p_kw);

% Effect size (eta-squared for ANOVA)
SS_between = tbl_anova{2,2};
SS_total = tbl_anova{4,2};
eta_squared = SS_between / SS_total;
fprintf('Effect size (η²): %.3f\n\n', eta_squared);

%%% 5. Post-hoc Pairwise Comparisons
if p_kw < alpha  % Use non-parametric post-hoc if omnibus significant
    fprintf('=== POST-HOC ANALYSIS (Non-parametric) ===\n');
    
    % Dunn's test (multiple comparisons after Kruskal-Wallis)
    [c, m, h, gnames] = multcompare(stats_kw, 'Display', 'off', 'CType', 'dunn-sidak');
    
    comparisons = {'Baseline vs Cooled', 'Baseline vs Recovery', 'Cooled vs Recovery'};
    
    for i = 1:size(c,1)
        fprintf('%s:\n', comparisons{i});
        fprintf('  Mean rank difference: %.2f\n', c(i,4));
        fprintf('  95%% CI: [%.2f, %.2f]\n', c(i,3), c(i,5));
        fprintf('  p-value: %.4f %s\n\n', c(i,6), ...
                ternary(c(i,6) < alpha, '(significant)', '(n.s.)'));
    end
end
%%% Box Plot Visualization
figure('Position', [100, 100, 800, 600]);

violinplot(all_data, group_labels);
title('Reaction Time Changes (Violin Plot)');
ylabel('% Change from Baseline Mean');
xlabel('Experimental Condition');



% Example p-values from your tests (replace with your actual results)
pvals = c(:,6); % [Baseline vs Cooled, Baseline vs Recovery, Cooled vs Recovery]


% Find y-limits for annotation placement
ylim_vals = ylim;
y_base = ylim_vals(2) + 0.05 * range(ylim_vals);

% Function to get significance stars

% Annotate pairwise comparisons
height = y_base;
step = 0.1 * range(ylim_vals);

% Baseline vs Cooled
line([1 1 2 2], [height height+step height+step height], 'Color', 'k', 'LineWidth', 1.5);
text(1.5, height+step, get_significance_stars(pvals(1)), 'HorizontalAlignment', 'center', 'FontSize', 14);

% Baseline vs Recovery
height = height + 1.5*step;
line([1 1 3 3], [height height+step height+step height], 'Color', 'k', 'LineWidth', 1.5);
text(2, height+step, get_significance_stars(pvals(2)), 'HorizontalAlignment', 'center', 'FontSize', 14);

% Cooled vs Recovery
height = height + 1.5*step;
line([2 2 3 3], [height height+step height+step height], 'Color', 'k', 'LineWidth', 1.5);
text(2.5, height+step, get_significance_stars(pvals(3)), 'HorizontalAlignment', 'center', 'FontSize', 14);

% Optionally, annotate the overall Kruskal-Wallis p-value
text(2, ylim_vals(1) + 0.95*range(ylim_vals), ...
    sprintf('Kruskal-Wallis p = %.3f', p_kw), ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');

hold off;



%% Statistical Summary
fprintf('Normalization Results:\n');
fprintf('Baseline Mean RT: %.3f ms\n', baseline_mean);
fprintf('\nMean Changes:\n');
fprintf('Baseline: %.2f%% ± %.2f%%\n', mean(baseline_norm), std(baseline_norm));
fprintf('Cooled: %.2f%% ± %.2f%%\n', mean(cooled_norm), std(cooled_norm));
fprintf('Recovery: %.2f%% ± %.2f%%\n', mean(recovery_norm), std(recovery_norm));

%% Plot out animal hit times and recovery period
load("\\10.165.57.13\Sutter_backup\Hammad\Ephys\LeverTask\Data_for_Figures\WavesCooling\M1CoolingWavesCombined.mat")


%% FUNCTIONS
function plotDat(dat)
baseline_data = dat(:,1);
eOPN_data = dat(:,2);
figure; % Create a new figure window

% Define colors for the points
baseline_color = [166 14 90]/255; % red
eOPN_color = [40 153 196]/255;   % cool blue

% Plot lines connecting paired points first (light gray, thin)
hold on; % Keep the plot active for multiple elements
for i = 1:numel(baseline_data)
    plot([1, 2], [baseline_data(i), eOPN_data(i)], 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5); % Light gray line
end

% Plot individual data points using scatter
% Baseline points (x=1)
scatter(ones(size(baseline_data)), baseline_data, 100, 'filled', ...
        'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', baseline_color);

% eOPN points (x=2)
scatter(2 * ones(size(eOPN_data)), eOPN_data, 100, 'filled', ...
        'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', eOPN_color);

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
xticklabels({'Baseline', 'Cooled'}); % Set x-axis labels
ax.XAxis.Color = [0.3 0.3 0.3]; % Darker color for x-axis labels
ax.YAxis.Color = [0.3 0.3 0.3]; % Darker color for y-axis labels
xlabel(''); % No overall x-axis label needed

% Y-axis label (match the image more closely)
ylabel('SI Index', 'FontSize', 16); 

ylim([0.0 1]); % Adjust limits to ensure 0 is visible

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
%% Helper function for effect size interpretation
function interp = effect_size_interpretation(d)
    if d < 0.2
        interp = '(negligible)';
    elseif d < 0.5
        interp = '(small)';
    elseif d < 0.8
        interp = '(medium)';
    else
        interp = '(large)';
    end
end

%% Helper function for ternary operator
function result = ternary(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end

function star_str = get_significance_stars(p_val)
    if p_val < 0.001
        star_str = '***';
    elseif p_val < 0.01
        star_str = '**';
    elseif p_val < 0.05
        star_str = '*';
    else
        star_str = 'n.s.';
    end
end