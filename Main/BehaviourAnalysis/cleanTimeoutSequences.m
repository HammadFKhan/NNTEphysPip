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
