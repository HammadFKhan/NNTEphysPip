
% Example usage:
sessionNums = [1]; % Example session numbers corresponding to days 1, 6, 12
sessionDays = [1]; % Corresponding training days for plotting

% Calculate metrics for each session
metrics = analyzeSequenceLearning(seqData, sessionNums);

% Plot metrics across days
plotLearningMetrics(metrics, sessionDays);
%%



% Function to analyze sequence learning metrics across sessions
function [metrics] = analyzeSequenceLearning(seqData, sessionNums)
    % Initialize structure to store metrics for each session
    metrics = struct();
    
    % Loop through selected sessions
    for i = 1:length(sessionNums)
        sessionNum = sessionNums(i);
        pullCounts = seqData.session(sessionNum).cleanedpullCounts;
        
        % Time axis (assuming time in ms, centered at reward)
        [numTrials, numTimePoints] = size(pullCounts);
        timeAxis = linspace(-1500, 1500, numTimePoints);
        timeStep = timeAxis(2) - timeAxis(1); % ms per time bin
        
        % Initialize arrays for metrics
        lengths = [];
        durations = [];
        isi = [];
        pressRates = [];
        
        % Loop through trials to extract sequences
        for trial = 1:numTrials
            % Find all non-zero entries (lever pulls)
            pullIdx = find(diff(pullCounts(trial, :)) > 0);
            pullTimes = timeAxis(pullIdx);
            pullValues = pullCounts(trial, pullIdx);
            
            % Consider only pulls before reward (negative times)
            preRewardIdx = pullTimes < 0;
            prePullTimes = pullTimes(preRewardIdx);
            prePullValues = pullValues(preRewardIdx);
            
            if ~isempty(prePullTimes)
                % Sequence length = number of pulls before reward
                seqLength = length(prePullTimes);
                lengths = [lengths; seqLength];
                
                % Sequence duration = time from first to last pull
                if seqLength > 1
                    seqDuration = abs(prePullTimes(end) - prePullTimes(1)) / 1000; % Convert to seconds
                    durations = [durations; seqDuration];
                    
                    % Within-sequence press rate (presses per minute)
                    pressRate = (seqLength - 1) / seqDuration * 60; % Convert to per minute
                    pressRates = [pressRates; pressRate];
                end
                
                % Calculate ISI (time between last press of previous trial and first press of current trial)
                if trial > 1
                    prevTrialPulls = find(pullCounts(trial-1, :) > 0);
                    if ~isempty(prevTrialPulls)
                        prevLastPullTime = timeAxis(prevTrialPulls(end));
                        currFirstPullTime = prePullTimes(1);
                        
                        % Convert to seconds and ensure positive value
                        trialISI = abs(currFirstPullTime - prevLastPullTime) / 1000;
                        isi = [isi; trialISI];
                    end
                end
            end
        end
        
        % Calculate metrics and store raw values
        metrics(i).session = sessionNum;

        % Store raw values arrays
        metrics(i).rawLengths = lengths;
        metrics(i).rawDurations = durations;
        metrics(i).rawISI = isi;
        metrics(i).rawPressRates = pressRates;

        % Calculate summary statistics
        metrics(i).meanLength = mean(lengths);
        metrics(i).cvLength = std(lengths) / mean(lengths);
        metrics(i).semLength = std(lengths) / sqrt(length(lengths));

        metrics(i).meanDuration = mean(durations);
        metrics(i).cvDuration = std(durations) / mean(durations);
        metrics(i).semDuration = std(durations) / sqrt(length(durations));

        metrics(i).meanISI = mean(isi);
        metrics(i).cvISI = std(isi) / mean(isi);
        metrics(i).semISI = std(isi) / sqrt(length(isi));

        metrics(i).meanPressRate = mean(pressRates);
        metrics(i).cvPressRate = std(pressRates) / mean(pressRates);
        metrics(i).semPressRate = std(pressRates) / sqrt(length(pressRates));

    end
end

% Function to plot learning metrics across sessions
function plotLearningMetrics(metrics, sessionDays)
    % Extract metric values
    meanLength = [metrics.meanLength];
    semLength = [metrics.semLength];
    cvLength = [metrics.cvLength];
    
    meanDuration = [metrics.meanDuration];
    semDuration = [metrics.semDuration];
    cvDuration = [metrics.cvDuration];
    
    meanISI = [metrics.meanISI];
    semISI = [metrics.semISI];
    cvISI = [metrics.cvISI];
    
    meanPressRate = [metrics.meanPressRate];
    semPressRate = [metrics.semPressRate];
    cvPressRate = [metrics.cvPressRate];
    
    % Create figure with 2 rows, 4 columns of subplots
    figure('Position', [100, 100, 1000, 500], 'Color', 'w');
    
    % Row 1: Means
    subplot(2,4,1);
    errorbar(sessionDays, meanLength, semLength, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    xlabel('Days of training');
    ylabel('Length of sequence (presses)');
    title('Sequence Length');
    
    subplot(2,4,2);
    errorbar(sessionDays, meanDuration, semDuration, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    xlabel('Days of training');
    ylabel('Duration of sequence (s)');
    title('Sequence Duration');
    
    subplot(2,4,3);
    errorbar(sessionDays, meanISI, semISI, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    xlabel('Days of training');
    ylabel('ISI (s)');
    title('Inter-Sequence Interval');
    
    subplot(2,4,4);
    errorbar(sessionDays, meanPressRate, semPressRate, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    xlabel('Days of training');
    ylabel('Within-sequence press rate (presses/min)');
    title('Press Rate');
    
    % Row 2: CVs
    subplot(2,4,5);
    errorbar(sessionDays, cvLength, semLength./meanLength, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    xlabel('Days of training');
    ylabel('Length of sequence (CV)');
    
    subplot(2,4,6);
    errorbar(sessionDays, cvDuration, semDuration./meanDuration, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    xlabel('Days of training');
    ylabel('Duration of sequence (CV)');
    
    subplot(2,4,7);
    errorbar(sessionDays, cvISI, semISI./meanISI, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    xlabel('Days of training');
    ylabel('ISI (CV)');
    
    subplot(2,4,8);
    errorbar(sessionDays, cvPressRate, semPressRate./meanPressRate, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    xlabel('Days of training');
    ylabel('Within-sequence press rate (CV)');
    
    % Add overall title
    sgtitle('Sequence Learning Metrics Across Training Days', 'FontSize', 16);
end
