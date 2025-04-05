%% Fit LM regression lines in figure 2
load('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\Trajectory\PooledDataNeuralSpeed.mat')
% Fitting time location of peak speed between cue and MI to reaction time
neuralTrajectoryTs = 20; % in ms
mdl = fitlm(neuralTrajectoryTs*PooledData.hit.speed.cueMIPeakSpeedTime,PooledData.RT)
figure,plot(mdl);
xlabel('Time of peak neural trajectory speed from cue');ylabel('Reaction Time');
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off

figure,scatter(neuralTrajectoryTs*PooledData.hit.speed.cueMIPeakSpeedTime,PooledData.RT,'filled','k');
xlabel('Time of peak neural trajectory speed from cue');ylabel('Reaction Time');
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
%%
% Fitting time location of peak speed between cue and MI to reaction time
mdl = fitlm(neuralTrajectoryTs*PooledData.hit.speed.cueMIPeakSpeed,PooledData.RT)
figure,plot(mdl);
xlabel('Peak neural trajectory speed from cue');ylabel('Reaction Time');
[ ~,outlierIndex] = rmoutliers(PooledData.RT);
PooledData.hit.speed.cueMIPeakSpeedRMOutlier = neuralTrajectoryTs*PooledData.hit.speed.cueMIPeakSpeed;
PooledData.hit.speed.cueMIPeakSpeedRMOutlier(outlierIndex) = [];
PooledData.RTRMOutlier = PooledData.RT; PooledData.RTRMOutlier(outlierIndex)=[];
mdl = fitlm(PooledData.hit.speed.cueMIPeakSpeedRMOutlier,PooledData.RTRMOutlier)
figure,plot(mdl);
xlabel('Peak neural trajectory speed from cue');ylabel('Reaction Time');
%%
load('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\Trajectory\PooledDataNeuralStability.mat')

figure,
scatter(stabilityPeak,reactionTime,'filled','k')
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
mdl = fitlm(stabilityPeak,reactionTime)
figure,plot(mdl)
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
ylim([0 1.5])