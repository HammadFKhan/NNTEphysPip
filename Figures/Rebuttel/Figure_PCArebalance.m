clear
redo = 1;
if redo==1
    files = dir(fullfile('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\M1_GSP','*.mat'));
    M1neuralDynamics = struct();
    for fileNum = 1:length(files)
        disp(['File number: ' num2str(fileNum)])
        load(fullfile(files(fileNum).folder,files(fileNum).name))
        try
            rebalanceTrialsforPCA = 0;
            [Spikes_unbalanced,M1neuralDynamics(fileNum).neuralDynamics_unbalanced] = getPCA(Spikes,rebalanceTrialsforPCA);
            rebalanceTrialsforPCA = 1;
            [Spikes_rebalanced,M1neuralDynamics(fileNum).neuralDynamics_rebalanced] = getPCA(Spikes,rebalanceTrialsforPCA);
            % Plot it all out
            M1neuralDynamics(fileNum).hitmiss_ub = Spikes_unbalanced.GPFA.resultHitMiss.kern.estParams.pcExplained*100;
            M1neuralDynamics(fileNum).hitfa_ub = Spikes_unbalanced.GPFA.resultMIHitFA.kern.estParams.pcExplained*100;

            M1neuralDynamics(fileNum).hitmiss_b = Spikes_rebalanced.GPFA.resultHitMiss.kern.estParams.pcExplained*100;
            M1neuralDynamics(fileNum).hitfa_b = Spikes_rebalanced.GPFA.resultMIHitFA.kern.estParams.pcExplained*100;
        catch ME
            disp(ME)
            disp('Error getting neural dynamics')
            continue
        end
        close all
    end
end

%% Make balanced PCA spikes
fileNum = 7;
hitmiss_ub = M1neuralDynamics(fileNum).hitmiss_ub;
hitmiss_b = M1neuralDynamics(fileNum).hitmiss_b;
hitfa_ub = M1neuralDynamics(fileNum).hitfa_ub;
hitfa_b = M1neuralDynamics(fileNum).hitfa_b;
figure,subplot(121),plotPCAexplained(hitmiss_ub),title('Unbalanced hitmiss')
subplot(122),plotPCAexplained(hitmiss_b),title('balanced hitmiss')

figure,subplot(121),plotPCAexplained(hitfa_ub),title('Unbalanced hitFA')
subplot(122),plotPCAexplained(hitfa_b),title('balanced hitFA')

%% stats across sessions
k = 5;  % or 5, 10
for n = 1:length(M1neuralDynamics)
    cum_ub = sum(M1neuralDynamics(n).hitmiss_ub(1:k), 2);
    cum_b  = sum(M1neuralDynamics(n).hitmiss_b(1:k), 2);

    [~, p_cum(n)] = ttest(cum_ub, cum_b);
    meanDiff_cum= mean(cum_ub - cum_b);
    meanDifftot(n,1) = meanDiff_cum;
    cum_ub = sum(M1neuralDynamics(n).hitfa_ub(1:k), 2);
    cum_b  = sum(M1neuralDynamics(n).hitfa_b(1:k), 2);
    meanDiff_cum= mean(cum_ub - cum_b);
    meanDifftot(n,2) = meanDiff_cum;
    [~, p_cum(n)] = ttest(cum_ub, cum_b);
    meanDiff_cum= mean(cum_ub - cum_b);
end
%%
figure,hold on,plotNiceBars(meanDifftot)
%%
fprintf('\nHitMiss Cumulative variance up to PC %d (unbalanced - balanced)\n', k);
fprintf('MeanDiff = %.2f%%, p = %.3g\n', ...
        meanDiff_cum, p_cum);

cum_ub = sum(hitfa_ub(1:k), 2);
cum_b  = sum(hitfa_b(1:k), 2);

[~, p_cum] = ttest(cum_ub, cum_b);
meanDiff_cum = mean(cum_ub - cum_b);
fprintf('\nHitFA Cumulative variance up to PC %d (unbalanced - balanced)\n', k);
fprintf('MeanDiff = %.2f%%, p = %.3g\n', ...
        meanDiff_cum, p_cum);
%% Take the larget balaance diff and plot it out the trajectory space
fileNum = 11;
neuralDynamics = M1neuralDynamics(fileNum).neuralDynamics_rebalanced;
x = squeeze(mean(neuralDynamics.hit.X(1,:,:),3));
y = squeeze(mean(neuralDynamics.hit.X(2,:,:),3));
z = squeeze(mean(neuralDynamics.hit.X(3,:,:),3));
figure,plot3(x,y,z),hold on
scatter3(x(75),y(75),z(75),25,'r','filled');
scatter3(x(88),y(88),z(88),25,'r','filled');
x = squeeze(mean(neuralDynamics.miss.X(1,:,:),3));
y = squeeze(mean(neuralDynamics.miss.X(2,:,:),3));
z = squeeze(mean(neuralDynamics.miss.X(3,:,:),3));
plot3(x,y,z,'m')
scatter3(x(75),y(75),z(75),25,'r','filled');
view(30,10)
set(gca,'tickdir','out')
axis square

x = squeeze(mean(neuralDynamics.MIhit.X(1,:,:),3));
y = squeeze(mean(neuralDynamics.MIhit.X(2,:,:),3));
z = squeeze(mean(neuralDynamics.MIhit.X(3,:,:),3));
figure,plot3(x,y,z),hold on
scatter3(x(75),y(75),z(75),25,'r','filled');
scatter3(x(88),y(88),z(88),25,'r','filled');
x = squeeze(mean(neuralDynamics.MIFA.X(1,:,:),3));
y = squeeze(mean(neuralDynamics.MIFA.X(2,:,:),3));
z = squeeze(mean(neuralDynamics.MIFA.X(3,:,:),3));
plot3(x,y,z,'m')
scatter3(x(75),y(75),z(75),25,'r','filled');
view(30,10)
set(gca,'tickdir','out')
axis square
%% Calaculate and plot neural speed across unbalanced datasets
dynamics = M1neuralDynamics;
dimension = 1;
speedTot = [];
for n = 1:length(dynamics)
speed_data = dynamics(n).neuralDynamics_unbalanced.hit.speed;
rt(n) = 12;
speedTot{n} = squeeze(speed_data.speed(dimension,2:end,:));
end
speedTot = horzcat(speedTot{:});

% figure;
% plot(speedTot);
% hold on;
% xline(75, '--r', 'Cue');
% xline(rt(n), '--g', 'Movement Start');
% xlabel('Time');
% ylabel('Average Speed');
% title(['Speed Over Time for Dimension ' num2str(dimension)]);
% legend('Speed', 'Cue', 'Movement Start');

colors = [0 0.4470 0.7410;0.75 0.75 0.75;190/255 30/255 45/255];
time = -1499:20:1500;
figure;
plot(time(2:end),mean(speedTot,2),'color',colors(1,:),'linewidth',2),hold on
plot(time(2:end),mean(speedTot,2)+std(speedTot,[],2)/sqrt(size(speedTot,2)),'color',colors(1,:),'linewidth',2)
plot(time(2:end),mean(speedTot,2)-std(speedTot,[],2)/sqrt(size(speedTot,2)),'color',colors(1,:),'linewidth',2)
hold on;

speedTot = [];
for n = 1:length(dynamics)
speed_data = dynamics(n).neuralDynamics_unbalanced.miss.speed;
rt(n) = 12 ;
speedTot{n} = squeeze(speed_data.speed(dimension,2:end,:));
end
rt(rt<10) = [];
speedTot = horzcat(speedTot{:});

plot(time(2:end),mean(speedTot,2),'color',colors(2,:),'linewidth',2),hold on
plot(time(2:end),mean(speedTot,2)+std(speedTot,[],2)/sqrt(size(speedTot,2)),'color',colors(2,:),'linewidth',2)
plot(time(2:end),mean(speedTot,2)-std(speedTot,[],2)/sqrt(size(speedTot,2)),'color',colors(2,:),'linewidth',2)

xline((75-75)*20, '--r', 'Cue');
xline((mean(rt)-75)*20, '--g', 'Movement Start');
xlabel('Time (s)');
ylabel('Average Speed');
box off,set(gca,'tickdir','out','fontsize',14),axis square,xlim([-500 1500]),ylim([0.01 0.1])

speedTot = [];
for n = 1:length(dynamics)
speed_data = dynamics(n).neuralDynamics_unbalanced.MIhit.speed;
rt(n) = 12  ;
speedTot{n} = squeeze(speed_data.speed(dimension,2:end,:));
end
speedTot = horzcat(speedTot{:});


colors = [0 0.4470 0.7410;0.75 0.75 0.75;190/255 30/255 45/255];
time = -1499:20:1500;
figure;
plot(time(2:end),mean(speedTot,2),'color',colors(1,:),'linewidth',2),hold on
plot(time(2:end),mean(speedTot,2)+std(speedTot,[],2)/sqrt(size(speedTot,2)),'color',colors(1,:),'linewidth',2)
plot(time(2:end),mean(speedTot,2)-std(speedTot,[],2)/sqrt(size(speedTot,2)),'color',colors(1,:),'linewidth',2)
hold on;
speedTot = [];

for n = 1:length(dynamics)
speed_data = dynamics(n).neuralDynamics_unbalanced.MIFA.speed;
rt(n) = 12  ;
speedTot{n} = squeeze(speed_data.speed(dimension,2:end,:));
end
speedTot = horzcat(speedTot{:});

plot(time(2:end),mean(speedTot,2),'color',colors(3,:),'linewidth',2),hold on
plot(time(2:end),mean(speedTot,2)+std(speedTot,[],2)/sqrt(size(speedTot,2)),'color',colors(3,:),'linewidth',2)
plot(time(2:end),mean(speedTot,2)-std(speedTot,[],2)/sqrt(size(speedTot,2)),'color',colors(3,:),'linewidth',2)

xline((75-75)*20, '--r', 'Movement Start');
xline((80-mean(rt))*20, '--g', 'Cue');
xlabel('Time (s)');
ylabel('Average Speed');
box off,set(gca,'tickdir','out','fontsize',14),axis square,xlim([-500 1500]),ylim([0.015 0.05])

%% Balance dataset
dynamics = M1neuralDynamics;
dimension = 1;
speedTot = [];
for n = 1:length(dynamics)
speed_data = dynamics(n).neuralDynamics_rebalanced.hit.speed;
rt(n) = 12;
speedTot{n} = squeeze(speed_data.speed(dimension,2:end,:));
end
speedTot = horzcat(speedTot{:});

% figure;
% plot(speedTot);
% hold on;
% xline(75, '--r', 'Cue');
% xline(rt(n), '--g', 'Movement Start');
% xlabel('Time');
% ylabel('Average Speed');
% title(['Speed Over Time for Dimension ' num2str(dimension)]);
% legend('Speed', 'Cue', 'Movement Start');

colors = [0 0.4470 0.7410;0.75 0.75 0.75;190/255 30/255 45/255];
time = -1499:20:1500;
figure;
plot(time(2:end),mean(speedTot,2),'color',colors(1,:),'linewidth',2),hold on
plot(time(2:end),mean(speedTot,2)+std(speedTot,[],2)/sqrt(size(speedTot,2)),'color',colors(1,:),'linewidth',2)
plot(time(2:end),mean(speedTot,2)-std(speedTot,[],2)/sqrt(size(speedTot,2)),'color',colors(1,:),'linewidth',2)
hold on;

speedTot = [];
for n = 1:length(dynamics)
speed_data = dynamics(n).neuralDynamics_rebalanced.miss.speed;
rt(n) = 12 ;
speedTot{n} = squeeze(speed_data.speed(dimension,2:end,:));
end
rt(rt<10) = [];
speedTot = horzcat(speedTot{:});

plot(time(2:end),mean(speedTot,2),'color',colors(2,:),'linewidth',2),hold on
plot(time(2:end),mean(speedTot,2)+std(speedTot,[],2)/sqrt(size(speedTot,2)),'color',colors(2,:),'linewidth',2)
plot(time(2:end),mean(speedTot,2)-std(speedTot,[],2)/sqrt(size(speedTot,2)),'color',colors(2,:),'linewidth',2)

xline((75-75)*20, '--r', 'Cue');
xline((mean(rt)-75)*20, '--g', 'Movement Start');
xlabel('Time (s)');
ylabel('Average Speed');
box off,set(gca,'tickdir','out','fontsize',14),axis square,xlim([-500 1500]),ylim([0.01 0.1])

speedTot = [];
for n = 1:length(dynamics)
speed_data = dynamics(n).neuralDynamics_rebalanced.MIhit.speed;
rt(n) = 12  ;
speedTot{n} = squeeze(speed_data.speed(dimension,2:end,:));
end
speedTot = horzcat(speedTot{:});


colors = [0 0.4470 0.7410;0.75 0.75 0.75;190/255 30/255 45/255];
time = -1499:20:1500;
figure;
plot(time(2:end),mean(speedTot,2),'color',colors(1,:),'linewidth',2),hold on
plot(time(2:end),mean(speedTot,2)+std(speedTot,[],2)/sqrt(size(speedTot,2)),'color',colors(1,:),'linewidth',2)
plot(time(2:end),mean(speedTot,2)-std(speedTot,[],2)/sqrt(size(speedTot,2)),'color',colors(1,:),'linewidth',2)
hold on;
speedTot = [];

for n = 1:length(dynamics)
speed_data = dynamics(n).neuralDynamics_rebalanced.MIFA.speed;
rt(n) = 12  ;
speedTot{n} = squeeze(speed_data.speed(dimension,2:end,:));
end
speedTot = horzcat(speedTot{:});

plot(time(2:end),mean(speedTot,2),'color',colors(3,:),'linewidth',2),hold on
plot(time(2:end),mean(speedTot,2)+std(speedTot,[],2)/sqrt(size(speedTot,2)),'color',colors(3,:),'linewidth',2)
plot(time(2:end),mean(speedTot,2)-std(speedTot,[],2)/sqrt(size(speedTot,2)),'color',colors(3,:),'linewidth',2)

xline((75-75)*20, '--r', 'Movement Start');
xline((80-mean(rt))*20, '--g', 'Cue');
xlabel('Time (s)');
ylabel('Average Speed');
box off,set(gca,'tickdir','out','fontsize',14),axis square,xlim([-500 1500]),ylim([0.015 0.05])

%% LOCAL FUNCTION

function [Spikes,neuralDynamics] = getPCA(Spikes,rebalanceTrialsforPCA)
Spikes = makeSpikeGPFA(Spikes,rebalanceTrialsforPCA);
Spikes.GPFA.HitMiss.dat = [Spikes.GPFA.hit.dat,Spikes.GPFA.miss.dat];
for n = length(Spikes.GPFA.hit.dat)+1:length(Spikes.GPFA.HitMiss.dat) %fix trials
    Spikes.GPFA.HitMiss.dat(n).trialId = n;
end
Spikes.GPFA.MIHitFA.dat = [Spikes.GPFA.MIHit.dat,Spikes.GPFA.MIFA.dat];
for n = length(Spikes.GPFA.MIHit.dat)+1:length(Spikes.GPFA.MIHitFA.dat) %fix trials
    Spikes.GPFA.MIHitFA.dat(n).trialId = n;
end
addpath(genpath('C:\Users\khan332\Documents\GitHub\NeuralTraj'));
addpath(genpath('mat_results'));
if exist('mat_results','dir'),rmdir('mat_results','s'),end
[Spikes.GPFA.resultHitMiss,Spikes.GPFA.seqTrainHitMiss] = gpfaAnalysis(Spikes.GPFA.HitMiss.dat,5); %Run index
[Spikes.GPFA.resultMIHitFA,Spikes.GPFA.seqTrainMIHitFA] = gpfaAnalysis(Spikes.GPFA.MIHitFA.dat,6); %Run index
neuralDynamics = neuralTrajBalance(Spikes);
end


function plotPCAexplained(pcaExplained)
% pcExplained: [nPC x 1] percent variance explained
nShow   = 6;                               % number of PCs to display
idx     = 1:nShow;
expVals = pcaExplained(idx);
cumVals = cumsum(pcaExplained);             % cumulative (%) over all PCs
cumVals = cumVals(idx);

hold on;

% 1) gray bars for individual PCs
bh = bar(idx, expVals, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');

% 2) teal cumulative line with black dots
ph = plot(idx, cumVals, '-o', ...
    'Color', [0 0.5 0.5], ...              % teal-ish line
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.5);

% 3) cosmetics
xlim([0.5 nShow+0.5]);
ylim([0 cumVals(end)+5]);
xlabel('PC');
ylabel('Variance Explained (%)');

set(gca, 'Box', 'off', ...
         'TickDir', 'out', ...
         'Layer', 'top','fontsize',10);
axis square
end

function neuralDynamics = neuralTrajBalance(Spikes)
if isfield(Spikes.GPFA,'seqTrainHitMiss')
X = arrayfun(@(x) vertcat(x.xorth),Spikes.GPFA.seqTrainHitMiss,'UniformOutput',false);
X = horzcat(X{:});
neuralTrajHitMiss = reshape(X,size(X,1),Spikes.GPFA.seqTrainHitMiss(1).T,[]);
end

if isfield(Spikes.GPFA,'seqTrainMIHitFA')
X = arrayfun(@(x) vertcat(x.xorth),Spikes.GPFA.seqTrainMIHitFA,'UniformOutput',false);
X = horzcat(X{:});
neuralTrajMIHitFA = reshape(X,size(X,1),Spikes.GPFA.seqTrainMIHitFA(1).T,[]);
end
dimNum = 15; %number of dimensions to take
if exist('neuralTrajHitMiss','var')
    X = neuralTrajHitMiss;
    hittrials = 1:length(Spikes.GPFA.hit.dat);
    misstrials = length(Spikes.GPFA.hit.dat)+1:size(X,3);

    [neuralDynamics.hit.r,neuralDynamics.hit.s,neuralDynamics.hit.stab,neuralDynamics.hit.X] = getMeanTraj(X,hittrials,dimNum); %trajectory variable and predefined conditional trial indexes
    neuralDynamics.hit.speed = speedTraj(X,hittrials,6);

    try
        [neuralDynamics.miss.r,neuralDynamics.miss.s,neuralDynamics.miss.stab,neuralDynamics.miss.X] = getMeanTraj(X,misstrials,dimNum); %trajectory variable and predefined conditional trial indexes

        neuralDynamics.miss.speed = speedTraj(X,misstrials,6);

        % Calculate differences in trajectories r'c(t)/||r'c(t)||
        [neuralDynamics.neuralSimhitmiss,neuralDynamics.neuralDiffhitmiss,neuralDynamics.hit.rprime,neuralDynamics.miss.rprime] = neuralTrajDiff(neuralDynamics.hit.r',neuralDynamics.miss.r');
    catch
        disp('error on neural dynamic miss')
    end
end

if exist('neuralTrajMIHitFA','var')
X = neuralTrajMIHitFA;
MIHittrials = 1:length(Spikes.GPFA.MIHit.dat);
MIFAtrials = length(Spikes.GPFA.MIHit.dat)+1:size(X,3);

[neuralDynamics.MIhit.r,neuralDynamics.MIhit.s,neuralDynamics.MIhit.stab,neuralDynamics.MIhit.X] = getMeanTraj(X,MIHittrials,dimNum); %trajectory variable and predefined conditional trial indexes
neuralDynamics.MIhit.speed = speedTraj(X,MIHittrials,6);


[neuralDynamics.MIFA.r,neuralDynamics.MIFA.s,neuralDynamics.MIFA.stab,neuralDynamics.MIFA.X] = getMeanTraj(X,MIFAtrials,dimNum); %trajectory variable and predefined conditional trial indexes
neuralDynamics.MIFA.speed = speedTraj(X,MIFAtrials,6);

% Calculate differences in trajectories r'c(t)/||r'c(t)||
[neuralDynamics.neuralSimMI,neuralDynamics.neuralDiffMI,neuralDynamics.MIhit.rprime,neuralDynamics.MIFA.rprime] = neuralTrajDiff(neuralDynamics.MIhit.r',neuralDynamics.MIFA.r');
end

end

function speed_struct = speedTraj(X,trials,components)
% data: 3D array (neural dimensions x time x trials)
% dimension: the neural dimension to analyze

% Extract the specified dimension
dimensions = 1:components;
dim_data = squeeze(X(dimensions,:, trials));

% Calculate the difference between consecutive time points
diff_data = diff(dim_data, 1, 2); % first derivative across the second dimension (time)

% Compute speed (absolute value of the difference)
speed = abs(diff_data);

% Add a row of zeros at the beginning to match original time dimension
speed_data = padarray(speed, [0 1 0], 0, 'pre');

% Calculate average speed for each state across all trials
% avg_speed_pre_cue = mean(speed_data(:,pre_cue_idx, :), [2 3]);
% avg_speed_pre_movement = mean(speed_data(:,pre_movement_idx, :), [2 3]);
% avg_speed_reward = mean(speed_data(:,reward_idx, :), [2 3]);

speed_struct.speed = speed_data;
end

function plotNiceBars(totData)
% totData: n x 6
% [nRows, nCols] = size(totData);
% npPoints = 24;
% repData = nan(npPoints, nCols);   % final 3 x 6 (3 points per column)
% for c = 1:nCols
%     x = totData(:, c);
%     x = x(~isnan(x));          % optional: drop NaNs per column
% 
%     mu = mean(x);
%     sd = std(x)/sqrt(length(x)*10);
% 
%     % target locations: mean, mean - sd, mean + sd
%     targets = [mu, mu - sd, mu + sd];
% 
%     % find indices of actual data closest to targets
%     idx = zeros(1, numel(targets));
%     for k = 1:numel(targets)
%         [~, idx(k)] = min(abs(x - targets(k)));
%     end
%     idx = unique(idx, 'stable');   % enforce uniqueness
% 
%     % if fewer than 3 unique points, fill remaining with random samples
%     if numel(idx) < npPoints
%         remaining = setdiff(1:numel(x), idx);
%         extra = randsample(remaining, npPoints - numel(idx));
%         idx = [idx, extra];
%     elseif numel(idx) > npPoints
%         idx = idx(1:npPoints);
%     end
% 
%     repData(:, c) = x(idx);
% end

% totData = repData;
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
set(gca, 'XTick', 1:size(totData,2),...
    'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
ylabel('IPI (s)');

hold off;

%%% RUN STATS
[p, tbl, stats] = anova1(totData, [], 'off'); % columns as groups
results = multcompare(stats, 'Display', 'off') % Pairwise comparisons

disp(['ANOVA p-value: ', num2str(p)]);
alpha = 0.05; % significance level
sigPairs = results(results(:,6) < alpha, :); % rows where p < 0.05
hold on;

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