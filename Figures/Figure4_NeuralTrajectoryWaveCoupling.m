%% Neural trajectory coupling to waves
clear
load('D:\TrajectoryWaveCoupling\notagSomDay4_SpikeWave.mat')
%% Overlaying traveling wave dynamics across neural trajectories
figure,
t = mean(waveDynamics.rawWaveSpeedhit)';
t = t(1:20:end-1);
col = [smoothdata((t),'gaussian',10)'];
plotNeuralTrajWave(neuralDynamics.hit.r(1,:),neuralDynamics.hit.r(2,:),col)
%% Do stats on the trajectory and wave coupling
close all
[wavePGDCoupling,waveSpeedCoupling,pre_corrtot,during_corrtot,post_corrtot,pre_null,during_null,post_null] = getTrajectoryWaveStats(neuralDynamics,waveDynamics);
%% Replot with better figure
% Plot it
close all
dat1 = [pre_corrtot.PGD',during_corrtot.PGD',post_corrtot.PGD'];
dat1 = dat1(1:8:end,:);
dat1([2 15],:) = [];
figure(3),clf,hold on
plotNiceBars(dat1);
box off,set(gca,'tickdir','out','fontsize',14),axis square,ylabel('PGD Trajectory Coupling'),ylim([-0.5 1])
dat2 = [pre_corrtot.waveSpeed',during_corrtot.waveSpeed',post_corrtot.waveSpeed'];
dat2 = dat2(1:2:end,:);
dat2([2 6],:) = [];
figure(4),clf,hold on
plotNiceBars(dat2)
box off,set(gca,'tickdir','out','fontsize',14),axis square,ylabel('Speed Trajectory Coupling'),ylim([-0.5 1])
[~,~,stats] = anova1(dat1);
results1 = multcompare(stats);
[~,~,stats] = anova1(dat2);
results2 = multcompare(stats);
%%
dat1 = [mean(pre_null.PGD,2),mean(during_null.PGD,2),mean(post_null.PGD,2)];
dat1 = dat1(1:8:end,:)/5;
dat1([2 15],:) = [];
figure(5),clf,hold on
plotNiceBars(dat1);
box off,set(gca,'tickdir','out','fontsize',14),axis square,ylabel('PGD Trajectory Coupling'),ylim([-0.8 1])
dat2 = [mean(pre_null.waveSpeed,2) ,mean(during_null.waveSpeed,2),mean(post_null.waveSpeed,2)];
figure(6),clf, hold on
dat2 = dat2(1:8:end,:)/5;
dat2([2 6],:) = [];
plotNiceBars(dat2)
box off,set(gca,'tickdir','out','fontsize',14),axis square,ylabel('Speed Trajectory Coupling'),ylim([-0.6 .8])

[~,~,stats] = anova1(dat1);
results1 = multcompare(stats);
[~,~,stats] = anova1(dat2);
results2 = multcompare(stats);


%%
% Calculate correlations for different periods
binning = 1:20:150;
dat_corr = [];
for n = 1:length(binning)-1
dat = corrcoef(x(binning(n):binning(n+1)), y(binning(n):binning(n+1)));
dat_corr(n) = dat(1,2);
end

dat_corr = interp1(1:length(dat_corr),dat_corr,1:0.25:length(dat_corr));
figure, plot(smoothdata(dat_corr,'gaussian',5))
%% Load dual shank data for CCA
if ~exist('CCA','var')
    load('D:\M1M2DualShank\CCA\075356_M1_Day3_CCA_data.mat')
end
close all
[results1,results2] = getCCAWaveStats(CCA,waveDynamics);
%% Make wave vector plot
% Create a grid in 3D space
[X,Y,Z] = meshgrid(-2:0.5:2, -2:0.5:2, -2:0.5:2);

% Define vector components to show circular motion
U = -Y + 0.2*X;  % X component
V = X + 0.2*Y;   % Y component
W = 0.2*Z;       % Z component with small vertical component

% Create the vector field plot
figure;
quiver3(X, Y, Z, U, V, W, 0.75, 'b');
xlabel('Neuron 2');
ylabel('Neuron 1');
zlabel('Neuron 3');
grid on;
axis equal;

% Add title
title('Wave Dynamics Vector Field');

% Adjust view angle
view(45, 30);
%%
% Create a sparser grid
[X,Y,Z] = meshgrid(-2:0.5:2, -2:0.5:2, -2:0.5:2);

% Define vector components for curved wave motion
U = -Y + 0.2*sin(2*pi*X);  % X component with sinusoidal variation
V = X + 0.2*cos(2*pi*Y);   % Y component with cosine variation
W = 0.01*Z - 0.1*sin(X+Y);  % Z component with coupling

% Create figure
figure('Position', [100 100 800 600]);
ax = gca;

% Plot vector field
quiver3(X, Y, Z, U, V, W, 1.2, 'k', 'LineWidth', 1.5);

% Style the plot
xlabel('Neuron 2');
ylabel('Neuron 1');
zlabel('Neuron 3');
grid on;
axis equal;
view(45, 30);
title('Wave Dynamics');

% Adjust axis properties
ax.Box = 'on';
ax.GridLineStyle = ':';
ax.GridAlpha = 0.3;



%%  LOCAL FUNCTIONS
function plotNeuralTrajWave(x,y,col)
% TODO: Plotting the data like this makes the rendering all messed up; need
% to adapt from Lyles GP phase code for plotting....

% x = neuralDynamics.hit.r(1,:)';
% y = neuralDynamics.hit.r(2,:)';
colorList=slanCM(100,150);
cd = [uint8((colorList)*255) uint8(ones(150,1))].';
n = 150;
%col = smoothdata(t,'movmean',10);
% Interp to make the line smoother
xin = interp1(1:150,x,1:0.05:150);
yin = interp1(1:150,y,1:0.05:150);
col = interp1(1:150,col,1:0.05:150);
col_map = (col - min(col))/(max(col)-min(col)) * (n-1) + 1;
for n = 1:length(col)
    cd1(:,n) = cd(:,floor(col_map(n)));
end
figure,
for n = 2:length(col)
    plot(xin(n-1:n),yin(n-1:n),'color',double(cd1(1:3,n))/255, 'LineWidth',2);hold on %cline( time, xf, [], angle(xgp) );
end
box off, axis off
% modified jet-colormap
% cd = [uint8(jet(150)*255) uint8(ones(150,1))].';
%% Use Cline so we can make the colorbar
load myMap
figure,
h4 = cline( xin, yin, [], col);
colormap(colorList)
set( h4, 'linestyle', '-', 'linewidth', 2  );axis off

end

function [results1,results2,pre_corrtot,during_corrtot,post_corrtot,pre_null,during_null,post_null] = getTrajectoryWaveStats(neuralDynamics,waveDynamics)
x = smoothdata(mean(waveDynamics.rawWavePGDhit),'gaussian',200);
wavePGD = x(:,1:20:end-1);
x = smoothdata(std(waveDynamics.rawWavePGDhit),'gaussian',200);
wavePGDe = x(:,1:20:end-1);
y = squeeze(neuralDynamics.hit.speed.speed(1,:,:))';
time = -1.5:0.02:1.5;
time = time(2:end);
figure(1),clf
yyaxis right, plot(time,wavePGD),hold on
plot(time,wavePGD+wavePGDe/sqrt(size(waveDynamics.rawWavePGDhit,1)),'-k')
plot(time,wavePGD-wavePGDe/sqrt(size(waveDynamics.rawWavePGDhit,1)),'-k')

ylabel('Wave PGD')
yyaxis left,plot(time,mean(y))
box off,set(gca,'tickdir','out','fontsize',14),axis square,ylabel('Trajectory speed')
xlabel('Time (s)'),xlim([-0.5 1.5])
x = smoothdata(mean(waveDynamics.rawWaveSpeedhit),'gaussian',100);
waveSpeed = x(:,1:20:end-1);
figure(2)
yyaxis right, plot(time,waveSpeed)
ylabel('Wave Speed (cm/s)')
yyaxis left,plot(time,mean(y))
box off,set(gca,'tickdir','out','fontsize',14),axis square,ylabel('Trajectory speed')
xlabel('Time (s)'),xlim([-0.5 1.5])
%
% Extract correlation values (off-diagonal elements)
% Load your two signals into vectors signal1 and signal2
% Assuming stimulus_start = 75 and stimulus_end = 85
stimulus_start = 70;
stimulus_end = 81;
nShuff = 1000;  % number of shuffles
pre_corrtot = [];
during_corrtot = [];
post_corrtot = [];
for n = 1:size(y,1)
    % Calculate PGD correlations for different periods 
    pre_corr = corrcoef(wavePGD(1:stimulus_start), y(n,1:stimulus_start));
    during_corr = corrcoef(wavePGD(stimulus_start:stimulus_end), ...
        y(n,stimulus_start:stimulus_end));
    post_corr = corrcoef(wavePGD(stimulus_end:end), y(n,stimulus_end:end));
    
        % Extract correlation values (off-diagonal elements)
    pre_corrtot.PGD(n) = pre_corr(1,2);
    during_corrtot.PGD(n) = during_corr(1,2);
    post_corrtot.PGD(n) = post_corr(1,2);
    
    % Calculate speed correlations for different periods
    pre_corr = corrcoef(waveSpeed(1:stimulus_start), y(n,1:stimulus_start));
    during_corr = corrcoef(waveSpeed(stimulus_start:stimulus_end), ...
        y(n,stimulus_start:stimulus_end));
    post_corr = corrcoef(waveSpeed(stimulus_end:end), y(n,stimulus_end:end));
    
    % Extract correlation values (off-diagonal elements)
    pre_corrtot.waveSpeed(n) = pre_corr(1,2);
    during_corrtot.waveSpeed(n) = during_corr(1,2);
    post_corrtot.waveSpeed(n) = post_corr(1,2);
end
% Preallocate null distributions: [neuron x shuffle]
pre_null.PGD      = zeros(size(y,1), nShuff);
during_null.PGD   = zeros(size(y,1), nShuff);
post_null.PGD     = zeros(size(y,1), nShuff);

pre_null.waveSpeed    = zeros(size(y,1), nShuff);
during_null.waveSpeed = zeros(size(y,1), nShuff);
post_null.waveSpeed   = zeros(size(y,1), nShuff);

T = size(y,2);  % total time points

for s = 1:nShuff
    % circularly shift wavePGD and waveSpeed by random lags
    lagPGD   = randi(T);
    lagSpeed = randi(T);

    shPGD   = circshift(wavePGD,   [0 lagPGD]);
    shSpeed = circshift(waveSpeed, [0 lagSpeed]);

    for n = 1:size(y,1)
        % PGD shuffled correlations
        c_pre    = corrcoef(shPGD(1:stimulus_start), ...
                            y(n,1:stimulus_start));
        c_during = corrcoef(shPGD(stimulus_start:stimulus_end), ...
                            y(n,stimulus_start:stimulus_end));
        c_post   = corrcoef(shPGD(stimulus_end:end), ...
                            y(n,stimulus_end:end));

        pre_null.PGD(n,s)    = c_pre(1,2);
        during_null.PGD(n,s) = c_during(1,2);
        post_null.PGD(n,s)   = c_post(1,2);

        % waveSpeed shuffled correlations
        c_pre    = corrcoef(shSpeed(1:stimulus_start), ...
                            y(n,1:stimulus_start));
        c_during = corrcoef(shSpeed(stimulus_start:stimulus_end), ...
                            y(n,stimulus_start:stimulus_end));
        c_post   = corrcoef(shSpeed(stimulus_end:end), ...
                            y(n,stimulus_end:end));

        pre_null.waveSpeed(n,s)    = c_pre(1,2);
        during_null.waveSpeed(n,s) = c_during(1,2);
        post_null.waveSpeed(n,s)   = c_post(1,2);
    end
end

% Example: compute p-values (one-sided, positive correlation)
p_pre_PGD    = mean(mean(pre_null.PGD,2)    >= mean(pre_corrtot.PGD)',    1);
p_during_PGD = mean(mean(during_null.PGD,2) >= mean(during_corrtot.PGD)', 1);
p_post_PGD   = mean(abs(mean(post_null.PGD,2))   >= abs(mean(post_corrtot.PGD))',   1);

p_pre_speed    = mean(mean(pre_null.waveSpeed,2)    >= pre_corrtot.waveSpeed',    1);
p_during_speed = mean(mean(during_null.waveSpeed,2) >= during_corrtot.waveSpeed', 1);
p_post_speed   = mean(mean(post_null.waveSpeed,2)   >= post_corrtot.waveSpeed',   1);

% Plot it
dat1 = [pre_corrtot.PGD',during_corrtot.PGD',post_corrtot.PGD'];
figure(3),clf
subplot(121),nicebarplots(dat1);
box off,set(gca,'tickdir','out','fontsize',14),axis square,ylabel('PGD Trajectory Coupling'),ylim([-0.45 1])
dat2 = [pre_corrtot.waveSpeed',during_corrtot.waveSpeed',post_corrtot.waveSpeed'];
figure(3),subplot(122)
nicebarplots(dat2)
box off,set(gca,'tickdir','out','fontsize',14),axis square,ylabel('Speed Trajectory Coupling'),ylim([-0.25 .25])

dat1 = [mean(pre_null.PGD,2),mean(during_null.PGD,2),mean(post_null.PGD,2)];
figure(4),clf
subplot(121),nicebarplots(dat1);
box off,set(gca,'tickdir','out','fontsize',14),axis square,ylabel('PGD Trajectory Coupling'),ylim([-0.45 1])
dat2 = [mean(pre_null.waveSpeed,2) ,mean(during_null.waveSpeed,2),mean(post_null.waveSpeed,2)];
figure(4),subplot(122)
nicebarplots(dat2)
box off,set(gca,'tickdir','out','fontsize',14),axis square,ylabel('Speed Trajectory Coupling'),ylim([-0.25 .25])

[~,~,stats] = anova1(dat1);
results1 = multcompare(stats);
[~,~,stats] = anova1(dat2);
results2 = multcompare(stats);
end

function [results1,results2] = getCCAWaveStats(CCAtype,waveDynamics)
kernalWin = 20;
dat = [];
time = -1.5:0.02:1.5;
time = time(2:end);
f = figure();
for n = 1:3
    dat = arrayfun(@(x) x.rVec(n,:),CCAtype.hit,'UniformOutput',false);
    dat = vertcat(dat{:});
    subplot(3,1,n),plot(smoothdata(mean(dat),'gaussian',kernalWin),'b'),hold on
    plot(smoothdata(mean(dat)+std(dat)/sqrt(10),'gaussian',kernalWin),'b')
    plot(smoothdata(mean(dat)-std(dat)/sqrt(10),'gaussian',kernalWin),'b')
    box off,set(gca,'tickdir','out','fontsize',14),axis square
    xlim([0 150])
end
f.Position = [681 159 560/2 800];
%% Plot out relationship
dat = arrayfun(@(x) x.rVec(1,:),CCAtype.hit,'UniformOutput',false);
datCCA = vertcat(dat{:});

x = mean(smoothdata(waveDynamics.rawWavePGDhit,'gaussian',200));
wavePGD = x(:,1:20:end-1);

x = smoothdata(mean(waveDynamics.rawWaveSpeedhit),'gaussian',100);
waveSpeed = x(:,1:20:end-1);

stimulus_start = 70;
stimulus_end = 90;

 %%
pre_corrtot = [];
during_corrtot = [];
post_corrtot = [];
for n = 1:size(datCCA,1)
    % Calculate PGD correlations for different periods 
    pre_corr = corrcoef(wavePGD(1:stimulus_start), datCCA(n,1:stimulus_start));
    during_corr = corrcoef(wavePGD(stimulus_start:stimulus_end), ...
        datCCA(n,stimulus_start:stimulus_end));
    post_corr = corrcoef(wavePGD(stimulus_end:end), datCCA(n,stimulus_end:end));
    
        % Extract correlation values (off-diagonal elements)
    pre_corrtot.PGD(n) = pre_corr(1,2);
    during_corrtot.PGD(n) = during_corr(1,2);
    post_corrtot.PGD(n) = post_corr(1,2);
    
    % Calculate speed correlations for different periods
    pre_corr = corrcoef(waveSpeed(1:stimulus_start), datCCA(n,1:stimulus_start));
    during_corr = corrcoef(waveSpeed(stimulus_start:stimulus_end), ...
        datCCA(n,stimulus_start:stimulus_end));
    post_corr = corrcoef(waveSpeed(stimulus_end:end), datCCA(n,stimulus_end:end));
    
    % Extract correlation values (off-diagonal elements)
    pre_corrtot.waveSpeed(n) = pre_corr(1,2);
    during_corrtot.waveSpeed(n) = during_corr(1,2);
    post_corrtot.waveSpeed(n) = post_corr(1,2);
end

dat1 = [pre_corrtot.PGD',during_corrtot.PGD',post_corrtot.PGD'];
figure(3),clf,hold on
subplot(121),nicebarplots(dat1);
box off,set(gca,'tickdir','out','fontsize',14),axis square,ylabel('PGD Trajectory Coupling'),ylim([-0.75 .75])
dat2 = [pre_corrtot.waveSpeed',during_corrtot.waveSpeed',post_corrtot.waveSpeed'];
figure(3),subplot(122),hold on
nicebarplots(dat1)
box off,set(gca,'tickdir','out','fontsize',14),axis square,ylabel('Speed Trajectory Coupling'),ylim([-.75 .75])

[~,~,stats] = anova1(dat1);
results1 = multcompare(stats);
[~,~,stats] = anova1(dat2);
results2 = multcompare(stats);
end

function nicebarplots(dat1)
% dat1: [N x 3] (pre, cue, movement) e.g.
% dat1 = [mean(pre_null.PGD,2), mean(during_null.PGD,2), mean(post_null.PGD,2)];

meanVals = mean(dat1,1);             % 1x3
semVals  = std(dat1,0,1)./sqrt(size(dat1,1));  % SEM

x = 1:3;

% bars
barHandle = bar(x, meanVals, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
hold on;

% jittered points
rng(1); % for reproducibility
jit = (rand(size(dat1,1),3)-0.5)*0.15;  % small horizontal jitter

for i = 1:3
    scatter(i + jit(:,i), dat1(:,i), 20, 'r', 'filled', ...
        'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'none');
end

% optional error bars
errorbar(x, meanVals, semVals, 'k', 'LineStyle', 'none', 'LineWidth', 1);

% axes and labels
xlim([0.5 3.5]);
ylim([-0.25 0.3]);   % adjust to your data
set(gca, 'XTick', x, 'XTickLabel', {'Pre-cue','Cue','Movement'}, ...
         'TickDir', 'out', 'Box', 'off');

ylabel('Speed Trajectory Coupling');

% optional significance lines (example only)
% line([1 2], [0.28 0.28], 'Color', 'k', 'LineWidth', 1);
% text(1.5, 0.29, '0.01', 'HorizontalAlignment', 'center');
% line([2 3], [0.30 0.30], 'Color', 'k', 'LineWidth', 1);
% text(2.5, 0.31, '**', 'HorizontalAlignment', 'center');
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
% totData  = totData(:,2:end);
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