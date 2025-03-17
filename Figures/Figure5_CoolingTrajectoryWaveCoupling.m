%% Figure 5 Trajectory wave coupling during cooling

load('Y:\Hammad\Ephys\LeverTask\DualShank\07538DualShankCoolingSOMCre\Day2\M130ChGridsDualShanksCoolingWavesDay2.mat')
load('Y:\Hammad\Ephys\LeverTask\DualShank\07538DualShankCoolingSOMCre\Day2\Day2DualShankGrid_CoolingRecording1_240918_120830\UCLA_chanmap_64F2\CCA_data.mat')

load('Y:\Hammad\Ephys\LeverTask\DualShank\07538DualShankCoolingSOMCre\Day8\M130ChGridsCoolingWavesDay8.mat')
load('Y:\Hammad\Ephys\LeverTask\DualShank\07538DualShankCoolingSOMCre\Day8\Shank64FM1GridM1CoolingM2Recording_240924_213957\UCLA_chanmap_64F2\Spikes.mat')
%% Sort waves cooled
tempCutoff = -9;
WavesBaseline = Waves;
WavesCooled = Waves;
WavesBaseline.wavesHit(IntanBehaviour.hitTemp<tempCutoff) = [];
WavesCooled.wavesHit(IntanBehaviour.hitTemp>=tempCutoff) = [];
%%
[M1neuralDynamics,M1waveDynamicsbaseline] = neuralTrajAnalysis2(Spikes,WavesBaseline,IntanBehaviour);
[~,M1waveDynamicscooled] = neuralTrajAnalysis2(Spikes,WavesCooled,IntanBehaviour);
%%% Sort neural dynamics
dynamicsBaseline = M1neuralDynamics;
dynamicsCooled = dynamicsBaseline;
temperatureId = (IntanBehaviour.hitTemp<tempCutoff);

dynamicsBaseline.hit.speed.speed = dynamicsBaseline.hit.speed.speed(:,:,~temperatureId);
dynamicsCooled.hit.speed.speed = dynamicsCooled.hit.speed.speed(:,:,temperatureId);
%% Recalculate trajectory subspace as a function of cooling
X = arrayfun(@(x) vertcat(x.xorth),Spikes.GPFA.seqTrainHitMiss,'UniformOutput',false);
X = horzcat(X{:});
neuralTrajHitMiss = reshape(X,size(X,1),Spikes.GPFA.seqTrainHitMiss(1).T,[]);

hitTrials = 1:length(IntanBehaviour.cueHitTrace);
hitTrials = hitTrials(~temperatureId);
[dynamicsBaseline.hit.r,dynamicsBaseline.hit.s] = meanTraj(neuralTrajHitMiss,hitTrials,6);
hitTrials = 1:length(IntanBehaviour.cueHitTrace);
hitTrials = hitTrials(temperatureId);
[dynamicsCooled.hit.r,dynamicsCooled.hit.s] = meanTraj(neuralTrajHitMiss,hitTrials,6);

%%
figure,
plot3(dynamicsBaseline.hit.r(1,:),dynamicsBaseline.hit.r(2,:),dynamicsBaseline.hit.r(3,:)),hold on
plot3(dynamicsCooled.hit.r(1,:),dynamicsCooled.hit.r(2,:),dynamicsCooled.hit.r(3,:))
%%
figure,
plotNeuralTrajWave(dynamicsBaseline.hit.r(1,:),dynamicsBaseline.hit.r(2,:),dynamicsBaseline.hit.r(3,:),M1waveDynamicsbaseline)
plotNeuralTrajWave(dynamicsCooled.hit.r(1,:),dynamicsCooled.hit.r(2,:),dynamicsCooled.hit.r(3,:),M1waveDynamicscooled)
%%
figure,
plotNeuralTrajWaveSpeed(dynamicsBaseline.hit.r(1,:),dynamicsBaseline.hit.r(2,:),dynamicsBaseline.hit.r(3,:),M1waveDynamicsbaseline)
plotNeuralTrajWaveSpeed(dynamicsCooled.hit.r(1,:),dynamicsCooled.hit.r(2,:),dynamicsCooled.hit.r(3,:),M1waveDynamicscooled)
%%
close all
[baselineDat,baselineStats] = getTrajectoryWaveStats(dynamicsBaseline,M1waveDynamicsbaseline);
[coolingDat,coolingStats] = getTrajectoryWaveStats(dynamicsCooled,M1waveDynamicscooled);
%% Calculate coupling difference
close all
% Check across subset of trials
PGDDiff = [];
SpeedDiff = [];
for n = 1:3
    trialIdBase = randperm(size(baselineDat.PGDTrajectorySpeed,1));
    trialIdBase = trialIdBase(1:floor(length(trialIdBase)*0.7));
    trialIdCool = randperm(size(coolingDat.PGDTrajectorySpeed,1));
    trialIdCool = trialIdCool(1:floor(length(trialIdCool)*0.7));
    PGDDiff(n,:) = abs((mean(baselineDat.PGDTrajectorySpeed(trialIdBase,:))-mean(coolingDat.PGDTrajectorySpeed(trialIdCool,:))))/2;
    SpeedDiff(n,:) = abs((mean(baselineDat.SpeedTrajectorySpeed(trialIdBase,:))-mean(coolingDat.SpeedTrajectorySpeed(trialIdCool,:))))/2;
end
figure,
subplot(121),errorbar(1:3, mean(PGDDiff,1),std(PGDDiff)*3)
xlim([0.5 3.5]),ylim([-0.1 0.7])

axis square,box off,set(gca,'tickdir','out','fontsize',14),axis square,ylabel('PGD Coupling difference')
subplot(122),errorbar(1:3, mean(SpeedDiff,1),std(SpeedDiff)*3),
xlim([0.5 3.5])
axis square,box off,set(gca,'tickdir','out','fontsize',14),axis square,ylabel('Speed Coupling difference')
ylim([-0.1 0.7])

%% LOCAL FUNCTIONS
function plotNeuralTrajWave(x,y,z, waveDynamics)
% TODO: Plotting the data like this makes the rendering all messed up; need
% to adapt from Lyles GP phase code for plotting....

% x = neuralDynamics.hit.r(1,:)';
% y = neuralDynamics.hit.r(2,:)';
cd = [uint8((jet(150))*255) uint8(ones(150,1))].';
n = 150;
t = mean(waveDynamics.rawWavePGDhit)';
t = t(1:20:end-1);
col = [0 smoothdata(abs(diff(t)),'gaussian',10)'];
%col = smoothdata(t,'movmean',10);
% Interp to make the line smoother
time = 1:20:3000;
xin = interp1(1:150,x,1:0.05:150);
yin = interp1(1:150,y,1:0.05:150);
zin = interp1(1:150,z,1:0.05:150);
col = interp1(1:150,col,1:0.05:150);
time = interp1(1:150,time,1:0.05:150);
col_map = (col - min(col))/(max(col)-min(col)) * (n-1) + 1;
for n = 1:length(col)
    cd1(:,n) = cd(:,floor(col_map(n)));
end
for n = 2:length(col)
    plot3(xin(n-1:n),yin(n-1:n),zin(n-1:n),'color',double(cd1(1:3,n))/255, 'LineWidth',2);hold on %cline( time, xf, [], angle(xgp) );
end
box off, axis off,view([0,0])
id = 1500; % cue
scatter3(xin(id),yin(id),zin(id),20,'k','filled')
id = 1850; % movement
scatter3(xin(id),yin(id),zin(id),20,'k','filled')
view([-10 40 ])
axis on
%axis([-0.5157    1.0006   -0.4182    0.0937   -0.2000    0.3000])
% modified jet-colormap
% cd = [uint8(jet(150)*255) uint8(ones(150,1))].';
%% Use Cline so we can make the colorbar
% load myMap
% figure,
% h4 = cline( xin, yin, [], col);
% colormap((jet))
% set( h4, 'linestyle', '-', 'linewidth', 2  );axis off

end
function plotNeuralTrajWaveSpeed(x,y,z, waveDynamics)
% TODO: Plotting the data like this makes the rendering all messed up; need
% to adapt from Lyles GP phase code for plotting....

% x = neuralDynamics.hit.r(1,:)';
% y = neuralDynamics.hit.r(2,:)';
cd = [uint8((jet(150))*255) uint8(ones(150,1))].';
n = 150;
t = arrayfun(@(x) mean(x.speed), waveDynamics.rawWaveSpeedHit)';
t = interp1(1:30,t,0:0.01:30);
t = inpaint_nans(t);
t = smoothdata(t,'gaussian',50)';
col = t(1:20:end-1);
%col = smoothdata(t,'movmean',10);
% Interp to make the line smoother
time = 1:20:3000;
xin = interp1(1:150,x,1:0.05:150);
yin = interp1(1:150,y,1:0.05:150);
zin = interp1(1:150,z,1:0.05:150);
col = interp1(1:150,col,1:0.05:150);
time = interp1(1:150,time,1:0.05:150);
col_map = (col - min(col))/(max(col)-min(col)) * (n-1) + 1;
for n = 1:length(col)
    cd1(:,n) = cd(:,floor(col_map(n)));
end
for n = 2:length(col)
    plot3(xin(n-1:n),yin(n-1:n),zin(n-1:n),'color',double(cd1(1:3,n))/255, 'LineWidth',2);hold on %cline( time, xf, [], angle(xgp) );
end
box off, axis off,view([0,0])
id = 1500; % cue
scatter3(xin(id),yin(id),zin(id),20,'k','filled')
id = 1850; % movement
scatter3(xin(id),yin(id),zin(id),20,'k','filled')
view([-10 40 ])
axis on
%axis([-0.5157    1.0006   -0.4182    0.0937   -0.2000    0.3000])
% modified jet-colormap
% cd = [uint8(jet(150)*255) uint8(ones(150,1))].';
%% Use Cline so we can make the colorbar
% load myMap
% figure,
% h4 = cline( xin, yin, [], col);
% colormap((jet))
% set( h4, 'linestyle', '-', 'linewidth', 2  );axis off

end

function [dat,stats] = getTrajectoryWaveStats(neuralDynamics,waveDynamics)
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
x = arrayfun(@(x) mean(x.speed),waveDynamics.rawWaveSpeedHit);
x = interp1(1:length(x),x,0:0.01:length(x));
x = smoothdata(x,'gaussian',500);
waveSpeed = x(:,1:20:end-1);
waveSpeed=inpaint_nans(waveSpeed);
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
figure(),clf
subplot(121),customBarplot(dat1);
box off,set(gca,'tickdir','out','fontsize',14),axis square,ylabel('PGD Trajectory Coupling'),ylim([-0.45 1])
dat2 = [pre_corrtot.waveSpeed',during_corrtot.waveSpeed',post_corrtot.waveSpeed'];
subplot(122)
customBarplot(dat2)
box off,set(gca,'tickdir','out','fontsize',14),axis square,ylabel('Speed Trajectory Coupling'),ylim([-0.25 .25])

[~,~,stats] = anova1(dat1);
results1 = multcompare(stats);
[~,~,stats] = anova1(dat2);
results2 = multcompare(stats);

stats = struct();
stats.result1 = results1;
stats.result2 = results2;
dat.PGDTrajectorySpeed = dat1;
dat.SpeedTrajectorySpeed = dat2;
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

function [r,s] = meanTraj(X,trials,components)
matrix = X(1:components,:,trials);
r = squeeze(mean(matrix,3));
% Subtract the mean from each element
centered_matrix = matrix - r;

% Calculate the Euclidean distance for each row
s = squeeze(sqrt(sum(centered_matrix.^2, 2)));

end
