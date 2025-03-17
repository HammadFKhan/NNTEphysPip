if ~exist('M1Spikes','var')||~exist('M2Spikes','var')
    disp('Loading data...')
    load('Y:\Hammad\Ephys\LeverTask\DualShank\075356DualShank\Day3\Day3M2DualM1SingleRecording1_240727_180222\CCA_data')
end

%% Sort by RT
% here we have the mean RT response from the normal data. To look at the
% effect of RT more closely we can take the 25th percentile above and below
% the mean to see how they change
[M1neuralDynamics,M1waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);
X = arrayfun(@(x) vertcat(x.xorth),Spikes.GPFA.seqTrainHitMiss,'UniformOutput',false);
X = horzcat(X{:});
neuralTrajHitMiss = reshape(X,size(X,1),Spikes.GPFA.seqTrainHitMiss(1).T,[]);
rttot = arrayfun(@(x) x.reactionTime,IntanBehaviour.cueHitTrace);
components = 6;

rtID = quantile(rttot,0.75); % Calculate quartile cutoff
hitTrials = 1:length(IntanBehaviour.cueHitTrace);
hitTrialsSlow = hitTrials(rttot>=rtID);
rtID = quantile(rttot,0.25); % Calculate quartile cutoff
hitTrialsFast = hitTrials(rttot<=rtID);

[rRTslow,~] = meanTraj(neuralTrajHitMiss,hitTrialsSlow,components);
[rRTfast,~] = meanTraj(neuralTrajHitMiss,hitTrialsFast,components);

%%% Plot it out
figure(1),clf
plot3(neuralDynamics.hit.r(1,:),neuralDynamics.hit.r(2,:),neuralDynamics.hit.r(3,:),'linewidth',2,'color',[0 0 0]),hold on
scatter3(neuralDynamics.hit.r(1,1:15:150),neuralDynamics.hit.r(2,1:15:150),neuralDynamics.hit.r(3,1:15:150),'filled','MarkerFaceColor',[0 0 0])

plot3(rRTslow(1,:),rRTslow(2,:),rRTslow(3,:),'linewidth',2,'color',[0.5 0.5 1])
scatter3(rRTslow(1,1:15:150),rRTslow(2,1:15:150),rRTslow(3,1:15:150),'filled','MarkerFaceColor',[0.5 0.5 1])

plot3(rRTfast(1,:),rRTfast(2,:),rRTfast(3,:),'linewidth',2,'color',[0 0 1])
scatter3(rRTfast(1,1:15:150),rRTfast(2,1:15:150),rRTfast(3,1:15:150),'filled','MarkerFaceColor',[0 0 1])

view([-45 30])
axis square
grid on
%% Plot distance as a function of reaction time
figure(2),clf,hold on
dim = 1;
dat1 = mean(neuralDynamics.hit.speed.speed(dim,:,hitTrialsSlow),3);
dat2 = mean(neuralDynamics.hit.speed.speed(dim,:,hitTrialsFast),3);
dat3 = mean(neuralDynamics.hit.speed.speed(dim,:,:),3);
plot(dat1,'LineWidth',2),plot(dat2,'LineWidth',2),plot(dat3,'LineWidth',2),legend('Slow','Fast','All')
%%
figure(3),clf,hold on
dim = 1;
dat1 = mean(neuralDynamics.hit.s(dim,hitTrialsSlow),2);
dat2 = mean(neuralDynamics.hit.s(dim,hitTrialsFast),2);
dat3 = mean(neuralDynamics.hit.s(dim,:),2);
bar([dat1,dat2,dat3])
%% Plot lever traces based on reactiontime
leverAll = arrayfun(@(x) x.trace, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
leverAll = horzcat(leverAll{:})';
leverSlow = leverAll(hitTrialsSlow,:);
leverFast = leverAll(hitTrialsFast,:);
figure(3),clf,hold
plot(mean(leverAll)),plot(mean(leverSlow)),plot(mean(leverFast)),legend('all','slow','fast')
figure,plot(mean(leverAll))
%%
colors = [0 0.4470 0.7410;0.75 0.75 0.75;190/255 30/255 45/255];
time = -1499:20:1500;
figure,plot(time(2:end),M1neuralDynamics.neuralDiffhitmiss(:,1:2),'color',colors(1,:),'LineWidth',2),ylim([-0.02 0.04]), hold on
plot(time(2:end),M1neuralDynamics.neuralDiffMI(:,1:2),'color',colors(2,:),'LineWidth',2),ylim([-0.02 0.04])
box off,set(gca,'tickdir','out','fontsize',14),axis square,xlim([-500 1500])
ylabel('Sequence Distance'),xlabel('Time (s)')

colors = [0 0.4470 0.7410;0.75 0.75 0.75;190/255 30/255 45/255];
time = -1499:20:1500;
figure,plot(time(2:end),M2neuralDynamics.neuralDiffhitmiss(:,1:2),'color',colors(1,:),'LineWidth',2),ylim([-0.06 0.16]), hold on
plot(time(2:end),M2neuralDynamics.neuralDiffMI(:,1:2),'color',colors(2,:),'LineWidth',2),ylim([-0.06 0.16])
box off,set(gca,'tickdir','out','fontsize',14),axis square,xlim([-500 1500])
ylabel('Sequence Distance'),xlabel('Time (s)')
%% Plot neural trajectory speed in relation to RT
neuralDynamics = M2neuralDynamics;
figure,scatter(neuralDynamics.PQ,neuralDynamics.rawreactionTime,'k','filled')
xlim([0 15])
ylim([0 2])
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
xlabel('Neural trajectory speed'),ylabel('Reaction time (s)')
mdl = fitlm(neuralDynamics.PQ,neuralDynamics.rawreactionTime)
figure,plot(mdl)
xlim([0 15])
ylim([0 2])
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
xlabel('Neural trajectory speed'),ylabel('Reaction time (s)')
legend(num2str(mdl.Rsquared.Ordinary),num2str(mdl.Coefficients.pValue(2)))
title('')
%% Plot neural trajectory STD
neuralDynamics = M1neuralDynamics;
temp = nan(max([size(neuralDynamics.hit.s,2),size(neuralDynamics.miss.s,2),size(neuralDynamics.MIFA.s,2)]),3);
temp(1:size(neuralDynamics.hit.s,2),1) = neuralDynamics.hit.s(1,:);
temp(1:size(neuralDynamics.miss.s,2),2) = neuralDynamics.miss.s(1,:);
temp(1:size(neuralDynamics.MIFA.s,2),3) = neuralDynamics.MIFA.s(1,:);

figure,customBarplot(temp);
box off,set(gca,'tickdir','out','fontsize',14),axis square
ylabel('Trajectory Deviation')
[p,t,stats] = anova1(temp)
c = multcompare(stats)
%% Trajectory Speed analysis
speed_data = M1neuralDynamics.hit.speed;
dimension = 1;
avg_speed_pre_cue = squeeze(mean(speed_data.preCue(dimension,:,:), [2]));
avg_speed_pre_movement = squeeze(mean(speed_data.preMovement(dimension,:,:), [2]));
avg_speed_reward = squeeze(mean(speed_data.reward(dimension,:,:), [2]));

temp = [avg_speed_pre_cue,avg_speed_pre_movement,avg_speed_reward];

% Plot average speed for each state
figure;
customBarplot(temp);
ylabel('Neural Trajectory Speed');
%title(['Average Trajectory Speed by Behavioral State for Dimension ' num2str(dimension)]);
box off,set(gca,'tickdir','out','fontsize',14),axis square
ylim([0 0.08])
%%
% Plot speed over time, highlighting different states
figure;
plot(mean(squeeze(speed_data.speed(1,:,:)), 2));
hold on;
xline(75, '--r', 'Cue');
xline(75+13, '--g', 'Movement Start');
xlabel('Time');
ylabel('Average Speed');
title(['Speed Over Time for Dimension ' num2str(dimension)]);
legend('Speed', 'Cue', 'Movement Start');

% Perform statistical comparison (ANOVA)
[p, tbl, stats] = anova1(temp);

fprintf('ANOVA p-value: %f\n', p);

if p < 0.05
    c = multcompare(stats, 'Display', 'off');
    disp('Multiple Comparisons:');
    disp(c);
end
%% LOCAL FUNCTIONS
function [r,s] = meanTraj(X,trials,components)
matrix = X(1:components,:,trials);
r = squeeze(mean(matrix,3));
% Subtract the mean from each element
centered_matrix = matrix - r;

% Calculate the Euclidean distance for each row
s = squeeze(sqrt(sum(centered_matrix.^2, 2)));

end