
parameters.experiment = 'self'; % self - internally generated, cue - cue initiated
parameters.opto = 0; % 1 - opto ON , 0 - opto OFF
parameters.cool = 0; % No Cool 
parameters.windowBeforePull = 1.5; % in seconds
parameters.windowAfterPull = 1.5; % in seconds
parameters.windowBeforeCue = 1.5; % in seconds
parameters.windowAfterCue = 1.5; % in seconds
parameters.windowBeforeMI = 1.5; % in seconds 
parameters.windowAfterMI = 1.5; % in seconds 
parameters.Fs = 1000; % Eventual downsampled data
parameters.ts = 1/parameters.Fs;

fname = 'Y:\Hammad\Ephys\SeqProject\Behavior\Rbp4M1Seq_Only2025_04_25_10.29.PM.csv';
[Behaviour] = readLeverSq(parameters);
%%
allPulls = arrayfun(@(x) x.pullCount, Behaviour.hitTrace, 'UniformOutput', false);
allPulls = horzcat(allPulls{:})';
%%% Make seqData structure
seqData = struct();
for n = 1
seqData.session(n).pullCounts = [zeros(size(allPulls,1),1),allPulls];
[seqData.session(n).cleanedpullCounts, hasTimeout] = cleanTimeoutSequences(seqData.session(n).pullCounts);
end
figure,imagesc(-1500:1500,1:184,seqData.session(n).cleanedpullCounts),ylabel('Trials'),xlabel('Time from reward (ms)')

%%

figure('Color', 'w', 'Position', [100, 100, 600, 900]);
imagesc([-1.5 1.5], [1 size(allPulls,1)], allPulls);

% Set color axis for clarity
caxis([0 3]);

% Add colorbar with label
cb = colorbar;
ylabel(cb, 'Lever Pull Count', 'FontSize', 14);

% Axis labels and title
xlabel('Time from Reward (ms)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Trial Number', 'FontSize', 14, 'FontWeight', 'bold');
title('Lever Pull Sequence Alignment to Reward', 'FontSize', 16, 'FontWeight', 'bold');

% Improve axis ticks
set(gca, 'FontSize', 12, 'XTick', -1.5:.5:1.5, 'YDir', 'normal');

% Optional: Add subtle gridlines
set(gca, 'YGrid', 'on', 'GridLineStyle', ':', 'GridAlpha', 0.2);

box off;
hold off;
%%
% Assume pullCounts is [trials x time] and timeAxis is the time vector (e.g., -1500:bin:1500)
allPulls = seqData.session(n).cleanedpullCounts;
[numTrials, numBins] = size(allPulls);
timeAxis = linspace(-3500, 500, numBins); % adjust as needed

durations = nan(numTrials,1);

for t = 1:numTrials
    % Only consider pulls before reward (time < 0)
    pulls = find(allPulls(t,:) > 0 & allPulls(t,:)< 4 & timeAxis < 0);
    if ~isempty(pulls)
        durations(t) = timeAxis(pulls(end)) - timeAxis(pulls(1));
    end
end

% Sort by duration (ascending: fastest to slowest)
[~, sortIdx] = sort(durations, 'ascend', 'MissingPlacement','last'); % NaNs (no pulls) go to bottom

% Reorder pullCounts for plotting
sortedPullCounts = allPulls(sortIdx, :);

% Plot
figure;
imagesc(timeAxis, 1:numTrials, sortedPullCounts);
colormap('parula'); % or your preferred colormap
colorbar;
xlabel('Time from reward (ms)');
ylabel('Trials (fastest to slowest)');
title('Lever Pulls Sorted by Sequence Duration');
%%
function [cleanedPullCounts, hasTimeout, sequenceStartedBeforeWindow] = cleanTimeoutSequences(pullCounts)
    [numTrials, numTimePoints] = size(pullCounts);
    cleanedPullCounts = zeros(size(pullCounts));
    hasTimeout = false(numTrials, 1);
    sequenceStartedBeforeWindow = false(numTrials, 1);
    timeAxis = linspace(-1500, 1500, numTimePoints);
    for trial = 1:numTrials
        % Get all non-zero entries (actual pull events)
        nonZeroMask = pullCounts(trial, :) > 0;
        timeIndices = find(nonZeroMask);
        
        if ~isempty(timeIndices)
            % Get the sequence of pull counts
            sequence = pullCounts(trial, timeIndices);
            pullTimes = timeAxis(timeIndices);
            
            % Find pulls before reward
            preRewardIndices = timeIndices(pullTimes < 0);
            
            if ~isempty(preRewardIndices)
                % Check if sequence starts with a value > 1 (started before window)
                firstValue = pullCounts(trial, preRewardIndices(1));
                if firstValue > 1
                    % Sequence was already in progress when recording started
                    sequenceStartedBeforeWindow(trial) = true;
                    % Keep trial zeroed out (exclude from analysis)
                    continue;
                end
                
                % Check for timeouts (non-monotonic counting)
                if length(sequence) > 1
                    diffs = diff(sequence);
                    
                    % If any difference is negative, we have a timeout/reset
                    if any(diffs < 0)
                        hasTimeout(trial) = true;
                        
                        % Find the last reset position
                        resetPositions = find(diffs < 0);
                        lastReset = resetPositions(end);
                        
                        % Start index of the final successful sequence
                        startIdx = timeIndices(lastReset + 1);
                        
                        % Set everything before final sequence to zero
                        cleanedPullCounts(trial, startIdx:end) = pullCounts(trial, startIdx:end);
                    else
                        % Monotonic sequence, no timeout
                        cleanedPullCounts(trial, :) = pullCounts(trial, :);
                    end
                else
                    % Only one pull, no possibility of reset
                    cleanedPullCounts(trial, :) = pullCounts(trial, :);
                end
            else
                % No pulls before reward, copy as is
                cleanedPullCounts(trial, :) = pullCounts(trial, :);
            end
        end
    end
    
    % Calculate statistics about excluded trials
    numTimeouts = sum(hasTimeout);
    numStartedBeforeWindow = sum(sequenceStartedBeforeWindow);
    totalExcluded = numTimeouts + numStartedBeforeWindow;
    
    fprintf('Trial statistics:\n');
    fprintf('  - Total trials: %d\n', numTrials);
    fprintf('  - Trials with timeouts: %d (%.1f%%)\n', numTimeouts, numTimeouts/numTrials*100);
    fprintf('  - Trials starting before window: %d (%.1f%%)\n', numStartedBeforeWindow, numStartedBeforeWindow/numTrials*100);
    fprintf('  - Total excluded trials: %d (%.1f%%)\n', totalExcluded, totalExcluded/numTrials*100);
end


