parameters.experiment = 'cue'; % self - internally generated, cue - cue initiated
parameters.opto = 0; % 1 - opto ON , 0 - opto OFF
parameters.cool = 0; % No Cool 
parameters.windowBeforePull = 3.5; % in seconds
parameters.windowAfterPull = 0.5; % in seconds
parameters.windowBeforeCue = 0.5; % in seconds
parameters.windowAfterCue = 3.5; % in seconds
parameters.windowBeforeMI = 1.5; % in seconds 
parameters.windowAfterMI = 1.5; % in seconds 
parameters.Fs = 1000; % Eventual downsampled data
parameters.ts = 1/parameters.Fs;
parameters.rows = 64;
parameters.cols = 1;
BehaviourSq = struct();
BehaviourSingle = struct();
fname = "Y:\Hammad\Ephys\SeqProject\Behavior\COM32025_05_11_03.01.PM.csv";
[BehaviourSq] = readLeverSqTrials(parameters,[],fname);
[BehaviourSingle] = readLeverSingleTrials(parameters,[],fname);
%%
figure,
for n = 1:length(BehaviourSingle.hitTrace)
    plot(smoothdata(BehaviourSingle.hitTrace(n).rawtrace),'color',[0.5 0.5 0.5 0.25]),hold on
end

figure,
for n = 1:length(BehaviourSq.hitTrace)
    plot(smoothdata(BehaviourSq.hitTrace(n).rawtrace),'color',[0.5 0.5 0.5 0.25]),hold on
end
%%

allPulls = arrayfun(@(x) x.pullCount, BehaviourSingle.cueHitTrace, 'UniformOutput', false);

% Determine the correct size (number of rows) from the first array
correctNumRows = size(allPulls{1}, 1);

% Find which arrays have the correct number of rows
validIdx = cellfun(@(c) size(c,1) == correctNumRows, allPulls);

% Keep only valid arrays
validPulls = allPulls(validIdx);

% Horizontally concatenate and transpose as you did
allPulls = horzcat(validPulls{:})';

[cleanedpullCounts, hasTimeout] = cleanTimeoutSequences(allPulls);
%%
data = cleanedpullCounts;
figure('Color', 'w', 'Position', [100, 100, 600, 900]);
imagesc([-3500 500], [1 size(data,1)], data);

% Use a perceptually uniform colormap

% Set color axis for clarity
caxis([0 3]);

% Add colorbar with label
cb = colorbar;
ylabel(cb, 'Lever Pull Count', 'FontSize', 14);

% Axis labels and title
xlabel('Time from Reward (ms)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Trial Number', 'FontSize', 14, 'FontWeight', 'bold');
title('Lever Pull Sequence Alignment to Reward', 'FontSize', 16, 'FontWeight', 'bold');

% Add vertical dashed line at reward delivery (time zero)
hold on;
plot([0 0], ylim, 'r--', 'LineWidth', 2);
text(30, size(data,1)-5, 'Reward', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold', 'Rotation', 90);
% 
% hold on;
% plot([-150 -150], ylim, 'r--', 'LineWidth', 2);
% text(30, size(data,1)-5, 'Reward', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold', 'Rotation', 90);

% Improve axis ticks
set(gca, 'FontSize', 12, 'XTick', -2500:500:500, 'YDir', 'normal');

% Optional: Add subtle gridlines
set(gca, 'YGrid', 'on', 'GridLineStyle', ':', 'GridAlpha', 0.2);

box off;
hold off;
