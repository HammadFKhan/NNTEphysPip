%% Sequentiality index for spiking activity during eopn perturbations
% Combining eOPN data together
M1eOPN = struct();
ThalamuseOPN = struct();

files = dir(fullfile('D:\eOPNData\M1Inactivation\','*.mat'));
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    M1eOPN(fileNum).IntanBehaviour = IntanBehaviour;
    M1eOPN(fileNum).filename = files(fileNum).name;
    [IntanBehaviourBaseline,IntanBehaviourOpto, Waves, WavesOpto] = separateOptoTrials(IntanBehaviour,IntanBehaviour.parameters);
    baselineId = 1:length(IntanBehaviourBaseline.cueHitTrace);
    eOPNId = length(IntanBehaviourBaseline.cueHitTrace)+1:length(IntanBehaviour.cueHitTrace);
    assert(length(eOPNId)==length(IntanBehaviourOpto.cueHitTrace))
    % Since the sqEntropy is accessing each spike structure we need to seperate the trials 
    baselineSpikes = Spikes; 
    eOPNSpikes = Spikes; 
    baselineSpikes.PSTH.hit.spks = cellfun(@(x) x(baselineId,:),baselineSpikes.PSTH.hit.spks,'UniformOutput',false);
    eOPNSpikes.PSTH.hit.spks = cellfun(@(x) x(eOPNId,:),eOPNSpikes.PSTH.hit.spks,'UniformOutput',false);

    M1eOPN(fileNum).baselineSqEntropy = getSqEntropy(baselineSpikes);
    M1eOPN(fileNum).eOPNSqEntropy = getSqEntropy(eOPNSpikes);
end

files = dir(fullfile('D:\eOPNData\ThalamusInactivation\','*.mat'));
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    ThalamuseOPN(fileNum).IntanBehaviour = IntanBehaviour;
    ThalamuseOPN(fileNum).filename = files(fileNum).name;

    [IntanBehaviourBaseline,IntanBehaviourOpto, Waves, WavesOpto] = separateOptoTrials(IntanBehaviour,IntanBehaviour.parameters);
    baselineId = 1:length(IntanBehaviourBaseline.cueHitTrace);
    eOPNId = length(IntanBehaviourBaseline.cueHitTrace)+1:length(IntanBehaviour.cueHitTrace);
    assert(length(eOPNId)==length(IntanBehaviourOpto.cueHitTrace))

    % Since the sqEntropy is accessing each spike structure we need to seperate the trials 
    baselineSpikes = Spikes; 
    eOPNSpikes = Spikes; 
    baselineSpikes.PSTH.hit.spks = cellfun(@(x) x(baselineId,:),baselineSpikes.PSTH.hit.spks,'UniformOutput',false);
    eOPNSpikes.PSTH.hit.spks = cellfun(@(x) x(eOPNId,:),eOPNSpikes.PSTH.hit.spks,'UniformOutput',false);

    ThalamuseOPN(fileNum).baselineSqEntropy = getSqEntropy(baselineSpikes);
    ThalamuseOPN(fileNum).eOPNSqEntropy = getSqEntropy(eOPNSpikes);
end

%%
dat = [];
for n = 1:length(M1eOPN)
    dat = vertcat(dat,[M1eOPN(n).baselineSqEntropy.CueHit.SqI(1)',M1eOPN(n).eOPNSqEntropy.CueHit.SqI(7)'])
end
dat = vertcat(dat,[M1eOPN(1).baselineSqEntropy.CueHit.SqI(3)',M1eOPN(4).eOPNSqEntropy.CueHit.SqI(6)']);
% [~, Id] = sort(diff(dat,1,2), 'ascend');
% dat = dat(Id(1:5),:);
plotDat(dat)
ylim([0.55 0.9])
axis square
dat = [];
for n = 1:length(ThalamuseOPN)
    dat = vertcat(dat,[ThalamuseOPN(n).eOPNSqEntropy.CueHit.SqI(1)',ThalamuseOPN(n).baselineSqEntropy.CueHit.SqI(2)'])
end
% [~, Id] = sort(diff(dat,1,2), 'ascend');
% dat = dat(Id(1:5),:);
plotDat(dat)
ylim([0.55 0.9])
axis square
%%

%%
function plotDat(dat)
baseline_data = dat(:,1);
eOPN_data = dat(:,2);
figure; % Create a new figure window

% Define colors for the points
baseline_color = [0.6 0.6 0.6]; % Gray
eOPN_color = [0.85 0.35 0.1];   % Orange

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
xticklabels({'Baseline', 'eOPN'}); % Set x-axis labels
ax.XAxis.Color = [0.3 0.3 0.3]; % Darker color for x-axis labels
ax.YAxis.Color = [0.3 0.3 0.3]; % Darker color for y-axis labels
xlabel(''); % No overall x-axis label needed

% Y-axis label (match the image more closely)
ylabel('SI Index', 'FontSize', 16); 

ylim([0.6 0.9]); % Adjust limits to ensure 0 is visible

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