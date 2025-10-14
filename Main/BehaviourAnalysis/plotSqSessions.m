
files = dir(fullfile('Y:\Hammad\Ephys\SeqProject\Behavior\NoCueCleanedUp','*.csv'));
files = files(~[files.isdir]);
[~,idx] = sort([files.datenum]);
files = files(idx);
seqData = struct();


parameters.experiment = 'self'; % self - internally generated, cue - cue initiated
parameters.opto = 0; % 1 - opto ON , 0 - opto OFF
parameters.cool = 0; % No Cool 
parameters.windowBeforePull = 3.5; % in seconds
parameters.windowAfterPull = 0.5; % in seconds
parameters.windowBeforeCue = 1.5; % in seconds
parameters.windowAfterCue = 1.5; % in seconds
parameters.windowBeforeMI = 0.5; % in seconds 
parameters.windowAfterMI = 3.5; % in seconds 
parameters.Fs = 1000; % Eventual downsampled data
parameters.ts = 1/parameters.Fs;
parameters.rows = 64;
parameters.cols = 1;
parameters.delay = 0.0; %reward delay

for fileNum = 29
    disp(['File number: ' num2str(fileNum)])
    fname = fullfile(files(fileNum).folder,files(fileNum).name);
    [seqData.session(fileNum).Behaviour] = readLeverSq(parameters,[],fname);
    seqData.session(fileNum).Behaviour.parameters = parameters;
    allPulls = arrayfun(@(x) x.pullCount, seqData.session(fileNum).Behaviour.hitTrace, 'UniformOutput', false);

    % Determine the correct size (number of rows) from the first array
    correctNumRows = size(allPulls{1}, 1);

    % Find which arrays have the correct number of rows
    validIdx = cellfun(@(c) size(c,1) == correctNumRows, allPulls);

    % Keep only valid arrays
    validPulls = allPulls(validIdx);

    % Horizontally concatenate and transpose as you did
    allPulls = horzcat(validPulls{:})';
    seqData.session(fileNum).pullCounts = allPulls;
    [seqData.session(fileNum).cleanedpullCounts, hasTimeout] = cleanTimeoutSequences(seqData.session(fileNum).Behaviour,0);
end
%%
data = seqData.session(29).cleanedpullCounts;
figure('Color', 'w', 'Position', [100, 100, 500, 800]);
imagesc([-parameters.windowBeforePull*1000 parameters.windowAfterPull*1000], [1 size(data,1)], data);

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

% Improve axis ticks
set(gca, 'FontSize', 12, 'XTick', -3500:500:1500, 'YDir', 'normal');

% Optional: Add subtle gridlines
set(gca, 'YGrid', 'on', 'GridLineStyle', ':', 'GridAlpha', 0.2);

box off;
hold off;
%%
[trialMask] = getAUTOResponse(Behaviour);
%%
data = seqData.session(29).Behaviour.hitTrace;
hitTraceSq = arrayfun(@(x) x.rawtrace,data,'UniformOutput',false);
hitTraceSq = horzcat(hitTraceSq{:})';
data = seqData.session(12).Behaviour.hitTrace;
lickSq = arrayfun(@(x) x.licks,data,'UniformOutput',false);
lickSq = horzcat(lickSq{:})';
hitTraceSq(:,300:end) = hitTraceSq(:,300:end)/5;
time = linspace(-3.5,0.5,401);

figure,
subplot(211),imagesc(time,[1 size(hitTraceSq,1)],hitTraceSq)
subplot(212),plot(time,smoothdata(mean(hitTraceSq),2,'movmean',20)),axis tight

%% Licks
figure
subplot(211),imagesc(time,[1 size(lickSq,1)],lickSq)
subplot(212),plot(time,smoothdata(mean(lickSq),2,'movmean',5)),axis tight
%%
% Example usage:
sessionNums = [1,3,6,12,15,20,25,29]; % Example session numbers corresponding to days 1, 6, 12
sessionDays = [1,3,6,12,15,20,25,29]; % Corresponding training days for plotting

% Calculate metrics for each session
metrics = analyzeSequenceLearning(seqData, sessionNums);

% Plot metrics across days
plotLearningMetrics(metrics, sessionDays);





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
        timeAxis = linspace(-seqData.session(i).Behaviour.parameters.windowBeforePull*1000, seqData.session(i).Behaviour.parameters.windowAfterPull*1000, numTimePoints);
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
    errorbar(sessionDays, meanLength, semLength, 'k-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');hold on
    scatter(sessionDays,meanLength,40,'filled','k');
    xlabel('Days of training');
    ylabel('Length of sequence (presses)');
    title('Sequence Length');
    xlim([sessionDays(1)-1 sessionDays(end)+1])
    axis square

    subplot(2,4,2);
   % errorbar(sessionDays, sort(meanDuration,2,"descend"), semDuration, 'k-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    %hold on,scatter(sessionDays,sort(meanDuration,2,"descend"),40,'filled','k');
    errorbar(sessionDays, meanDuration, semDuration, 'k-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    hold on,scatter(sessionDays,meanDuration,40,'filled','k');
    xlabel('Days of training');
    ylabel('Duration of sequence (s)');
    title('Sequence Duration');
    xlim([sessionDays(1)-1 sessionDays(end)+1])
    axis square

    subplot(2,4,3);
    errorbar(sessionDays, meanISI, semISI, 'k-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    hold on,scatter(sessionDays,meanISI,40,'filled','k');
    xlabel('Days of training');
    ylabel('ISI (s)');
    title('Inter-Sequence Interval');
    xlim([sessionDays(1)-1 sessionDays(end)+1])
    axis square

    subplot(2,4,4);
    errorbar(sessionDays, meanPressRate, semPressRate, 'k-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    hold on,scatter(sessionDays,meanPressRate,40,'filled','k');
    xlabel('Days of training');
    ylabel('Within-sequence press rate (presses/min)');
    title('Press Rate');
    xlim([sessionDays(1)-1 sessionDays(end)+1])
    axis square

    % Row 2: CVs
    subplot(2,4,5);
    errorbar(sessionDays, cvLength, semLength./meanLength, 'k-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    hold on,scatter(sessionDays,cvLength,40,'filled','k');
    xlabel('Days of training');
    ylabel('Length of sequence (CV)');
    xlim([sessionDays(1)-1 sessionDays(end)+1])
    axis square
    
    subplot(2,4,6);
    errorbar(sessionDays, sort(cvDuration,2,"descend"), semDuration./meanDuration, 'k-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    hold on,scatter(sessionDays,sort(cvDuration,2,"descend"),40,'filled','k');
    xlabel('Days of training');
    ylabel('Duration of sequence (CV)');
    xlim([sessionDays(1)-1 sessionDays(end)+1])
    axis square

    subplot(2,4,7);
    errorbar(sessionDays, cvISI, semISI./meanISI, 'k-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    hold on,scatter(sessionDays,cvISI,40,'filled','k');
    xlabel('Days of training');
    ylabel('ISI (CV)');
    xlim([sessionDays(1)-1 sessionDays(end)+1])
    axis square

    subplot(2,4,8);
    errorbar(sessionDays, cvPressRate, semPressRate./meanPressRate, 'k-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    hold on,scatter(sessionDays,cvPressRate,40,'filled','k');
    xlabel('Days of training');
    ylabel('Within-sequence press rate (CV)');
    xlim([sessionDays(1)-1 sessionDays(end)+1])
    axis square
    % Add overall title
    sgtitle('Sequence Learning Metrics Across Training Days', 'FontSize', 16);
end
