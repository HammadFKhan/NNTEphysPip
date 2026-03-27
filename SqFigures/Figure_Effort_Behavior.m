%% Combine and produce lever Intan files based on overall animal behavior during force field perturbation task. 
% Directory to load in Intan Behavior data
files = dir(fullfile('D:\SQLever\Ephys\ForceField\','*.mat'));
IntanBehaviourPooled = struct();
for fileNum = 1:length(files)
    fName = fullfile(files(fileNum).folder,files(fileNum).name);
    disp(['Loading ' fName '...'])
    load(fName)
    IntanBehaviourPooled(fileNum).IntanBehaviour = IntanBehaviour;
    IntanBehaviourPooled(fileNum).fname = files(fileNum).name;
end
%% Plot out effort and lever traces
time = linspace(0,5,1001);
f = figure;
f.Position = [680 458 660 520];
subplot(2,2,1)
for i=1:length(IntanBehaviour.MIHitTrace)
    plot(time,smoothdata(IntanBehaviour.MIHitTrace(i).trace(1:5:end)),'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
end
box off,set(gca,'tickdir','out')
axis square,xlim([0 5])
subplot(2,2,3)
for i=1:length(IntanBehaviour.hitTrace)
    plot(time,(IntanBehaviour.hitTrace(i).effortTrace(1:5:end)),'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
end
box off,set(gca,'tickdir','out')
axis square,xlim([0 5])

subplot(2,2,2)
for i=1:length(IntanBehaviour.effortperturbTrace)
    plot(time,smoothdata(IntanBehaviour.effortperturbTrace(i).trace(1:5:end)),'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
end
box off,set(gca,'tickdir','out')
axis square,xlim([0 5])

subplot(2,2,4)
for i=1:length(IntanBehaviour.effortperturbTrace)
    plot(time,(IntanBehaviour.effortperturbTrace(i).effortTrace(1:5:end)),'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
end
box off,set(gca,'tickdir','out')
axis square,xlim([0 5])
%% Calculate inter-pull interval
% IntanBehaviour = IntanBehaviourPooled(1).IntanBehaviour;
pullIPI = vertcat(IntanBehaviour.hitTrace.pullCount);
pullIPI = diff(pullIPI,1,2)/1000;
figure;subplot(121), hold on;
cleanDat = pullIPI;
cleanDat = rmoutliers(cleanDat,'quartiles');
plotNiceBars(cleanDat),ylim([0 1.6])

% Make IPI for perturb pull counts
pullIPI = nan(length(IntanBehaviour.effortperturbTrace),IntanBehaviour.SqNum-1);
for n = 1:length(IntanBehaviour.effortperturbTrace)
    pullTrial = IntanBehaviour.effortperturbTrace(n).pullCount;
    % Condition if there is a breakthrough on the catch trials (ie. more
    % than 1 pull)
    pullTrial(pullTrial<IntanBehaviour.parameters.windowBeforeMI*IntanBehaviour.parameters.Fs) = [];
    if length(pullTrial)>1
        temp = diff(pullTrial)/1000;
        % Add condition where the third pull was not completed (Sq didn't
        % finish
        if length(temp)<2
            temp(2) = NaN;
        end
        if length(temp)>2
            temp = temp(1:2);
        end
        pullIPI(n,:) = temp;
    end
end

subplot(122), hold on;
plotNiceBars(pullIPI),ylim([0 1.6])
%% Calculate distance that the lever was pulled
% We consider the z=scored lever amp as the peak lever amplitude from the
% pull count index to the proceeding pull count. Assumption is based on the
% a threshold crossing to find the max peak. 

[normative_LA,catch_LA] = getLeverAmp(IntanBehaviourPooled(1).IntanBehaviour);
figure,subplot(121),hold on
plotNiceBars(normative_LA)
set(gca, 'XTick', 1:size(normative_LA,2), 'XTickLabel', {'1', '2', '3', 'Late'}, ...
    'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
ylabel('Lever Distance');
ylim([0 6]);
subplot(122),hold on
plotNiceBars(catch_LA)
set(gca, 'XTick', 1:size(catch_LA,2), 'XTickLabel', {'1', '2', '3', 'Late'}, ...
    'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
ylabel('Lever Distance');
ylim([0 6]);
%% Calculate probability of breakout trials
% We can calculate the brekout trials from the total number of pulls
% generated from the leverAmp response
% Breakout prob is the number of trials outside the first pull. Which is
% only the second column onward. But we can also calculated a weighted
% pullout if we want to treat the third column (ie. completed sequence) as
% more important
for n = 1:length(IntanBehaviourPooled)
    [normative_LA,catch_LA] = getLeverAmp(IntanBehaviourPooled(n).IntanBehaviour);
    pullOut = sum(~isnan(catch_LA));
    breakoutProb(n)= pullOut(3)/pullOut(1);
    weighted_breakoutprob(n) = sum(pullOut(2:end))/pullOut(1);
    % Calculate lever amplitude shift out
    normLA(n,:) = nanmean(normative_LA,1);
    norlmLAe(n,:) = nanstd(normative_LA,[],1);
    catchLA(n,:) = nanmean(catch_LA,1);
    catchLAe(n,:) = nanstd(catch_LA,[],1);
end
%%%
% Example bins (adjust as needed)
sessions = length(breakoutProb);
x = 1:sessions;
% Example error bars (replace with your actual errors if available)
err1 = 0; % or your SEM/CI
err2 = 0; % or your SEM/CI

figure;
% breakoutProb = sort(breakoutProb);
% weighted_breakoutprob = sort(weighted_breakoutprob);
% Top Panel
subplot(2,1,1); hold on;
fill([x fliplr(x)], [breakoutProb+err1 fliplr(breakoutProb-err1)], [0.8 0.6 0.8], 'EdgeColor','none', 'FaceAlpha',0.4);
plot(x, breakoutProb, '-o', 'Color', [0.7 0.3 0.7], 'MarkerFaceColor', [0.7 0.3 0.7], 'LineWidth',2);
set(gca, 'XTickLabel', [], 'XColor', 'none');
ylabel('Breakthrough trials');
xlim([0 length(breakoutProb)+1]); set(gca,'Box','off', 'FontSize',12,'tickdir','out');

% Bottom Panel
subplot(2,1,2); hold on;
fill([x fliplr(x)], [weighted_breakoutprob+err2 fliplr(weighted_breakoutprob-err2)], [0.5 0.2 0.7], 'EdgeColor','none', 'FaceAlpha',0.4);
plot(x, weighted_breakoutprob, '-o', 'Color', [0.5 0.2 0.7], 'MarkerFaceColor', [0.5 0.2 0.7], 'LineWidth',2);


xlabel('Sessions');
ylabel('Weighted Breakthrough trials');
xlim([0 length(breakoutProb)+1]); set(gca,'Box','off', 'FontSize',12,'tickdir','out');

% Example bins (adjust as needed)

% Example error bars (replace with your actual errors if available)
err1 = 0; % or your SEM/CI
err2 = 0; % or your SEM/CI

%%% Plot Lever Amp

sessions = size(normLA,1);
x = 1:sessions; % Session (or time, or bin) indices
nLines = size(normLA,2);
lineColors = lines(nLines); % Distinct colors for each line
norlmLAe = norlmLAe/sqrt(sessions);
catchLAe = catchLAe/sqrt(sessions);
catchLAe = fillmissing(catchLAe, 'linear');
catchLA = fillmissing(catchLA, 'linear');

baseColor1 = [0.7 0.3 0.7]; % base purple for top plot

baseColor2 = [0.5 0.2 0.7]; % base purple for bottom plot

figure;

% Top panel
subplot(2,1,1); hold on;
for k = 1:nLines
    shadeFactor = 0.6 + (k-1)*0.3; % e.g. 0.6, 0.4, 0.2
    lineColor = baseColor1 * shadeFactor; 
    % Shaded error
    fill([x fliplr(x)], ...
        [normLA(:,k)'+norlmLAe(:,k)' fliplr(normLA(:,k)'-norlmLAe(:,k)')], ...
        lineColor, 'FaceAlpha',0.2, 'EdgeColor','none');
    % Mean line with markers
    plot(x, normLA(:,k), '-o', 'Color', lineColor, 'MarkerFaceColor',lineColor, 'LineWidth',2);
end
ylabel('Lever amplitude');
set(gca, 'XTickLabel', [], 'XColor', 'none', 'Box','off', 'FontSize', 12,'tickdir','out'); % No x-axis line/labels
xlim([0 length(breakoutProb)+1]);
ylim([0 6])
% Bottom panel

subplot(2,1,2); hold on;
for k = 1:nLines
    shadeFactor = 0.6 + (k-1)*0.4; % e.g. 0.6, 0.4, 0.2
    lineColor = baseColor2 * shadeFactor; 
    fill([x fliplr(x)], ...
         [catchLA(:,k)'+catchLAe(:,k)' fliplr(catchLA(:,k)'-catchLAe(:,k)')], ...
         lineColor, 'FaceAlpha',0.2, 'EdgeColor','none');
    plot(x, catchLA(:,k), '-o', 'Color', lineColor, 'MarkerFaceColor', lineColor, 'LineWidth',2);
end
xlabel('Session');
ylabel('Lever amplitude');
set(gca, 'Box','off', 'FontSize', 12,'tickdir','out');
xlim([0 length(breakoutProb)+1]);
ylim([0 6])
%%
force = [7;25];
figure;hold on
scatter([1,2],force,25,'filled'),xlim([0 2])
plot([1,2], force, '-', ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 1); % semi-transparent gray
axis square
xlim([0.5 ,2.5])
ylim([0 30])
set(gca, 'XTick', 1:2, 'XTickLabel', {'Normative', 'Catch'}, ...
    'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
ylabel('Pull Weight (g)')
%% Build a subsequent trial analysis of the effort response
% Here we need to figure out the chronological order of the effort task vs
% the hit task....
% We can do this by sorting the index of the pull count to the LFP time for
% both effort and hit trials
periCatchTrajectorySessions = [];
periCatchSpeedSessions =[];
leverSpeedNorm = [];
leverSpeedCatch = [];
for nSess = 1:length(IntanBehaviourPooled)
    IntanBehaviour = IntanBehaviourPooled(nSess).IntanBehaviour;
    effortPullTime = nan(1,length(IntanBehaviour.effortperturbTrace));
    hitTrialPullTime = nan(1,length(IntanBehaviour.MIHitTrace));
    for n = 1:length(IntanBehaviour.hitTrace)
        trueTime = IntanBehaviour.MIHitTrace(n).LFPIndex;
        firstPull = IntanBehaviour.MIHitTrace(n).pullCount(1);
        hitTrialPullTime(n) = trueTime(firstPull);
    end
    for n = 1:length(IntanBehaviour.effortperturbTrace)
        trueTime = IntanBehaviour.effortperturbTrace(n).LFPIndex;
        firstPull = IntanBehaviour.effortperturbTrace(n).pullCount(1);
        effortPullTime(n) = trueTime(firstPull);
    end

    % Find the index for shared trials (ie. when a breakout trial was marked
    % out)
    % It also is required because of some bugs we have....
    BO_trials = intersect(hitTrialPullTime,effortPullTime);
    if ~isempty(BO_trials)
        for BO = 1:length(BO_trials)
            trialID = find(hitTrialPullTime==BO_trials(BO));
            IntanBehaviour.hitTrace(trialID).effortFlag = 1;
        end
        %assert(length(BO_trials)==length(vertcat(IntanBehaviour.effortperturbTrace.rewardFlag)))
    end

    allTrials = [hitTrialPullTime,effortPullTime];
    allTrials(2,:) = [ones(1,length(hitTrialPullTime)),zeros(1,length(effortPullTime))];
    allTrials(3,:) = [1:length(hitTrialPullTime),1:length(effortPullTime)];
    [~,c] = sort(allTrials(1,:)); %Sort chronologically
    allTrials = allTrials(:,c);
    %%% Trajectory correlation
    % We can calculate the lever trajectory correlation between normative catch
    % trials as described by Shadmehr and Mussa-ilvaldi. Importantly we want to
    % make a triggered trial window around the catch trial _=2 trials, which we
    % assume is the normative trials
    catchTrials = find(allTrials(2,:)==0);
    normTrials = find(allTrials(2,:)==1);
    periCatchTrajectory = nan(length(catchTrials),3);
    periCatchSpeed = nan(length(catchTrials),3);
    leverTrace = smoothdata([horzcat(IntanBehaviour.MIHitTrace.trace),horzcat(IntanBehaviour.effortperturbTrace.trace)],'gaussian',5)+min(horzcat(IntanBehaviour.MIHitTrace.trace),[],'all');
    leverTrace = (diff(leverTrace+10));
    %figure,plot(leverTrace)
    u = mean(leverTrace(:,1:length(IntanBehaviour.MIHitTrace)),2);
    leverTrace  = leverTrace(:,c);
    catchTrials(end) = [];
    %     leverSpeedCatch{nSess} = max(leverTrace(:,catchTrials));
    %     leverSpeedNorm{nSess} = max(leverTrace(:,catchTrials-1));
    %     leverSpeedNorm2{nSess} = max(leverTrace(:,catchTrials+1));
    % for loop below seems to work quite well
%     for n = -1:1
%         leverSpeedCatch{nSess,n+2} = max(leverTrace(:,catchTrials+n))';
%     end
    % But not here...
    win = -2:2;
    for n = 1:length(catchTrials)
        for trial = win
            idx = catchTrials(n) + trial;
            if idx < 1 || idx > size(leverTrace, 2)
                continue; % skip out of range indices
            end
            try
                v = leverTrace(:, idx);  % Use same index logic as working loop
                periCatchSpeed(n, trial + abs(min(win))+1) = max(v);
                CC = cov(v,u)/(std(v)*std(u));
                periCatchTrajectory(n,trial+abs(min(win))+1) = CC(1,2);
            catch
                disp('error computing periCatchSpeed...');
                continue;
            end
        end
    end
    periCatchTrajectorySessions{nSess} = periCatchTrajectory;
    periCatchSpeedSessions{nSess} =  periCatchSpeed;
end
%% Plot session average
dat = vertcat(periCatchTrajectorySessions{:});
figure;
plot(win, abs(dat)', 'Color', [0.7 0.7 0.7], 'LineWidth', 1); % All catch trials, gray
hold on
plot(win, nanmean(abs(dat),1), '-o', 'Color',[0.7 0.3 0.7], 'MarkerFaceColor',[0.7 0.3 0.7], 'LineWidth', 1.5); % Mean, thick black line
set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlabel('Trial from perturbation');
ylabel('CC');
xlim([-2.5 2.5]);
ylim([0 0.65]);
axis square
hold off
%% 
dat = vertcat(periCatchSpeedSessions{:});
figure;
plot(win, abs(dat)', 'Color', [0.7 0.7 0.7], 'LineWidth', 1); % All catch trials, gray
hold on
plot(win, nanmean(abs(dat),1), '-o', 'Color',[0.7 0.3 0.7], 'MarkerFaceColor',[0.7 0.3 0.7], 'LineWidth', 1.5); % Mean, thick black line
set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlabel('Trial from perturbation');
ylabel('Lever Speed');
xlim([-2.5 2.5]);
axis square
%% plot raw relationship of lever speed
dat = vertcat(leverSpeedCatch{:,1});
figure,customBarplot([vertcat(leverSpeedCatch{:,1}),vertcat(leverSpeedCatch{:,2}),vertcat(leverSpeedCatch{:,3})])
lm = fitlm(x, y);
% Display the model summary
disp(lm);

% Plot the data and the fitted regression line
figure;
scatter(x, y, 'filled');
hold on;
plot(lm);
xlabel('Mean leverSpeedNorm');
ylabel('Mean leverSpeedCatch');
title('Linear Regression Model Fit with fitlm');
hold off;

%%
function plotNiceBars(totData)
means = nanmean(totData);          % Bar heights
sems = nanstd(totData) ./ sqrt(size(totData,1));   % Error bar (standard error)
b = bar(means, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k'); % Gray bars with black edge

% Overlay error bars
errorbar(1:size(totData,2), means, sems, 'k', 'LineStyle', 'none', 'LineWidth', 1);

% Overlay individual jittered points
xjitter = randn(size(totData))*0.01; % Controls point jitter
for i = 1:size(totData,2)
    scatter(i + xjitter(:,min(i,2)), totData(:,i), 18, 'o', ...
        'MarkerEdgeColor', [0.25 0.25 0.25], ...
        'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4);
end

% Draw paired lines between columns 1 and 2
for j = 1:size(totData,1)
    if size(totData,2) >= 3  % If there are at least 3 columns
        xvals = [1 + xjitter(j,1), 2 + xjitter(j,2), 3 + xjitter(j,3)];
        yvals = [totData(j,1),    totData(j,2),    totData(j,3)];
        plot(xvals, yvals, '-', 'Color', [0.5 0.5 0.5 0.6], 'LineWidth', 1);
    else % Connect just columns 1 and 2
        xvals = [1 + xjitter(j,1), 2 + xjitter(j,2)];
        yvals = [totData(j,1),    totData(j,2)];
        plot(xvals, yvals, '-', 'Color', [0.5 0.5 0.5 0.6], 'LineWidth', 1);
    end
end

% Style similar to image
set(gca, 'XTick', 1:size(totData,2), 'XTickLabel', {'Second Pull', 'Third Pull', 'Polymer', 'Late'}, ...
    'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
ylabel('IPI (s)');
ylim([0 2]);

hold off;

%%% RUN STATS
[p, tbl, stats] = anova1(totData, [], 'off'); % columns as groups
results = multcompare(stats, 'Display', 'off'); % Pairwise comparisons

disp(['ANOVA p-value: ', num2str(p)]);
alpha = 0.05; % significance level
sigPairs = results(results(:,6) < alpha, :); % rows where p < 0.05
hold on;
ylims = ylim;

% vertical height offset for significance lines above bars
baseY = max(means + sems) * 1.05;  
offsetStep = max(means + sems) * 0.05; 
if all(results(:,6) >= 0.05) % No significant pairwise differences
    % Extract F statistic from ANOVA table
    Fstat = cell2mat(tbl(2,5)); % Assumes standard anova1 output tbl
    p_anova = p;
    % Place text on plot upper corner
    xPos = size(totData,2)/2;
    yPos = max(means + sems) * 2.4;
    text(xPos, yPos, sprintf('ANOVA F=%.2f, p=%.3f', Fstat, p_anova), ...
        'HorizontalAlignment', 'left', 'FontSize', 10);
    % Add pairwise stars or p-values as before (your existing code)
end

for i = 1:size(sigPairs,1)
    x1 = sigPairs(i,1);
    x2 = sigPairs(i,2);
    y = baseY + (i-1)*offsetStep;
    
    % Draw line connecting bars
    plot([x1 x1 x2 x2], [y y+offsetStep y+offsetStep y], 'k-', 'LineWidth', 1);
    
    % Add star above the line
    text(mean([x1 x2]), y + offsetStep*0.1, '*', 'HorizontalAlignment', 'center', ...
        'FontSize', 16, 'FontWeight', 'bold');
end
hold off;
end

function [normative_LA,catch_LA] = getLeverAmp(IntanBehaviour)
leverAmp = nan(length(IntanBehaviour.hitTrace),IntanBehaviour.SqNum);
leverTrace = smoothdata(horzcat(IntanBehaviour.hitTrace.trace))+min(horzcat(IntanBehaviour.hitTrace.trace),[],'all');
%leverTrace = (leverTrace-min(leverTrace,[],'all'))/(max(leverTrace,[],'all')-min(leverTrace,[],'all'));
leverTrace = zscore(leverTrace);
for n = 1:length(IntanBehaviour.hitTrace)
    for nn = 1:length(IntanBehaviour.hitTrace(n).pullCount)
        if nn == length(IntanBehaviour.hitTrace(n).pullCount)
            leverAmp(n,nn) =  max(leverTrace(IntanBehaviour.hitTrace(n).pullCount(nn):end,n));
        else
            leverAmp(n,nn) =  max(leverTrace(IntanBehaviour.hitTrace(n).pullCount(nn):IntanBehaviour.hitTrace(n).pullCount(nn+1),n));
        end
    end
end
normative_LA = leverAmp;

leverAmp = nan(length(IntanBehaviour.effortperturbTrace),IntanBehaviour.SqNum);
leverTrace = smoothdata(horzcat(IntanBehaviour.effortperturbTrace.trace))+min(horzcat(IntanBehaviour.effortperturbTrace.trace),[],'all');
%leverTrace = (leverTrace-min(leverTrace,[],'all'))/(max(leverTrace,[],'all')-min(leverTrace,[],'all'));
leverTrace = zscore(leverTrace);
for n = 1:length(IntanBehaviour.effortperturbTrace)
    pullTrial = IntanBehaviour.effortperturbTrace(n).pullCount;
    pullTrial(pullTrial<IntanBehaviour.parameters.windowBeforeMI*IntanBehaviour.parameters.Fs) = [];
    if length(pullTrial)>IntanBehaviour.SqNum
        pullTrial = pullTrial(1:IntanBehaviour.SqNum);
    end
    for nn = 1:length(pullTrial)
        if nn == length(pullTrial)
            leverAmp(n,nn) =  max(leverTrace(pullTrial(nn):end,n));
        else
            leverAmp(n,nn) =  max(leverTrace(pullTrial(nn):pullTrial(nn+1),n));
        end
    end
end
catch_LA = leverAmp;
end