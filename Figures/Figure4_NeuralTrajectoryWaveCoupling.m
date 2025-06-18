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
[wavePGDCoupling,waveSpeedCoupling] = getTrajectoryWaveStats(neuralDynamics,waveDynamics);
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

function [results1,results2] = getTrajectoryWaveStats(neuralDynamics,waveDynamics)
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
% Plot it
dat1 = [pre_corrtot.PGD',during_corrtot.PGD',post_corrtot.PGD'];
figure(3),clf
subplot(121),customBarplot(dat1);
box off,set(gca,'tickdir','out','fontsize',14),axis square,ylabel('PGD Trajectory Coupling'),ylim([-0.45 1])
dat2 = [pre_corrtot.waveSpeed',during_corrtot.waveSpeed',post_corrtot.waveSpeed'];
figure(3),subplot(122)
customBarplot(dat2)
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
figure(3),clf
subplot(121),customBarplot(dat1);
box off,set(gca,'tickdir','out','fontsize',14),axis square,ylabel('PGD Trajectory Coupling'),ylim([-0.75 .75])
dat2 = [pre_corrtot.waveSpeed',during_corrtot.waveSpeed',post_corrtot.waveSpeed'];
figure(3),subplot(122)
customBarplot(dat2)
box off,set(gca,'tickdir','out','fontsize',14),axis square,ylabel('Speed Trajectory Coupling'),ylim([-.75 .75])

[~,~,stats] = anova1(dat1);
results1 = multcompare(stats);
[~,~,stats] = anova1(dat2);
results2 = multcompare(stats);
end