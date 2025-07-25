function [cleanedPullCounts, pullIndices,Behaviour,hasTimeout, sequenceStartedBeforeWindow] = cleanTimeoutSequences(Behaviour,MIFlag)


if ~exist('MIFlag','var')
    MIFlag = 1;
    disp('Cleaning up MI flags')
    allPulls = arrayfun(@(x) x.pullCount, Behaviour.MIhitTrace, 'UniformOutput', false);
end
if MIFlag == 0
    disp('Cleaning up hit trace')
    allPulls = arrayfun(@(x) x.pullCount, Behaviour.hitTrace, 'UniformOutput', false);
end
% Determine the correct size (number of rows) from the first array
correctNumRows = size(allPulls{1}, 1);

% Find which arrays have the correct number of rows
validIdx = cellfun(@(c) size(c,1) == correctNumRows, allPulls);

% Keep only valid arrays
validPulls = allPulls(validIdx);

% Horizontally concatenate and transpose as you did
pullCounts = horzcat(validPulls{:})';
[numTrials, numTimePoints] = size(pullCounts);
cleanedPullCounts = zeros(size(pullCounts));
hasTimeout = false(numTrials, 1);
sequenceStartedBeforeWindow = false(numTrials, 1);
pullIndices = cell(numTrials,1);
% Check if we want to cleanup based on movement initiation or reward
% aligned data
if MIFlag == 1
    timeAxis = linspace(-Behaviour.parameters.windowBeforeMI, Behaviour.parameters.windowAfterMI  , numTimePoints);
    for trial = 1:numTrials
        % Get all non-zero entries (actual pull events)
        nonZeroMask = pullCounts(trial, :) > 0;
        timeIndices = find(nonZeroMask);

        if ~isempty(timeIndices)
            % Get the sequence of pull counts
            sequence = pullCounts(trial, timeIndices);
            pullTimes = timeAxis(timeIndices);
            % Find pulls before reward
            preRewardIndices = pullTimes < Behaviour.MIHitTrace(trial).rewardtime;
            timeIndices = timeIndices(preRewardIndices);
            sequence = sequence(preRewardIndices);
            if ~isempty(preRewardIndices)
                % Check if sequence starts with a value > 1 (started before window)
                firstValue = sequence(1);
                if firstValue > 1
                    % Sequence was already in progress when recording started
                    sequenceStartedBeforeWindow(trial) = true;
                    % Keep trial zeroed out (exclude from analysis)
                    continue;
                end

                % Check for timeouts (non-monotonic counting)
                if length(unique(sequence)) > 0
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
        cleanedPullCounts(trial,(Behaviour.MIHitTrace(trial).rewardtime*100+151):end) = zeros(1,length(cleanedPullCounts(trial,(Behaviour.MIHitTrace(trial).rewardtime*100+151):end)));
        Behaviour.MIHitTrace(trial).cleanedPullCounts = cleanedPullCounts(trial,:)';
    end
else
    timeAxis = linspace(-Behaviour.parameters.windowBeforePull, Behaviour.parameters.windowAfterPull  , numTimePoints);
    for trial = 1:numTrials
        % Get all non-zero entries (actual pull events)
        nonZeroMask = pullCounts(trial, :) > 0;
        timeIndices = find(nonZeroMask);

        if ~isempty(timeIndices)
            % Get the sequence of pull counts
            sequence = pullCounts(trial, timeIndices);
            pullTimes = timeAxis(timeIndices);
            % Find pulls before reward
            preRewardIndices = pullTimes < -Behaviour.parameters.delay  ; %reward time is at zero
            timeIndices = timeIndices(preRewardIndices);
            sequence = sequence(preRewardIndices);
            if ~isempty(preRewardIndices)
                % Check if sequence starts with a value > 1 (started before window)
                firstValue = sequence(1);
                if firstValue > 1
                    % Sequence was already in progress when recording started
                    sequenceStartedBeforeWindow(trial) = true;
                    % Keep trial zeroed out (exclude from analysis)
                    continue;
                end

                % Check for timeouts (non-monotonic counting)
                if length(unique(sequence)) > 0
                    diffs = diff(sequence);
                    pullStart = find(diffs==1)+1;
                    pullIndices{trial} = [timeIndices(1),timeIndices(pullStart)];
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
        cleanedPullCounts(trial,find(timeAxis==0):end) = zeros(1,length(cleanedPullCounts(trial,find(timeAxis==0):end)));
        Behaviour.hitTrace(trial).cleanedPullCounts = cleanedPullCounts(trial,:)';
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
