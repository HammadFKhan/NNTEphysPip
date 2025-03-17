clear
%fpath = 'F:\LeverTask\Ephys\Analysis\M2Spikes';
%fpath = 'F:\LeverTask\Ephys\Analysis\spksPooledwFA'
%fpath = 'D:\M1_GSP';
fpath = 'D:\M2SpikeData';
file = dir(fullfile(fpath,'*.mat'));
M1Dynamics = struct();

for fileNum = 1:length(file)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(file(fileNum).folder,file(fileNum).name))
    Spikes = makeSpikeGPFA(Spikes);
    Spikes.GPFA.HitMiss.dat = [Spikes.GPFA.hit.dat,Spikes.GPFA.miss.dat];
    for n = 1:IntanBehaviour.nCueHit%+1:length(Spikes.GPFA.BaselineOpto.dat) %fix trials
        Spikes.GPFA.HitMiss.dat(n).trialId = n;
    end
    Spikes.GPFA.MIHitFA.dat = [Spikes.GPFA.MIHit.dat,Spikes.GPFA.MIFA.dat];
    for n = length(IntanBehaviour.MIHitTrace)+1:length(Spikes.GPFA.MIHitFA.dat) %fix trials
        Spikes.GPFA.MIHitFA.dat(n).trialId = n;
    end
    %%%
    addpath(genpath('C:\Users\khan332\Documents\GitHub\NeuralTraj'));
    addpath(genpath('mat_results'));
    if exist('mat_results','dir'),rmdir('mat_results','s'),end
    [Spikes.GPFA.resultHit,Spikes.GPFA.seqTrainHit] = gpfaAnalysis(Spikes.GPFA.hit.dat,1); %Run index
    [Spikes.GPFA.resultMiss,Spikes.GPFA.seqTrainMiss] = gpfaAnalysis(Spikes.GPFA.miss.dat,2); %Run index
    [Spikes.GPFA.resultMIHit,Spikes.GPFA.seqTrainMIHit] = gpfaAnalysis(Spikes.GPFA.MIHit.dat,3); %Run index
    [Spikes.GPFA.resultMIFA,Spikes.GPFA.seqTrainMIFA] = gpfaAnalysis(Spikes.GPFA.MIFA.dat,4); %Run index
    [Spikes.GPFA.resultHitMiss,Spikes.GPFA.seqTrainHitMiss] = gpfaAnalysis(Spikes.GPFA.HitMiss.dat,5); %Run index
    [Spikes.GPFA.resultMIHitFA,Spikes.GPFA.seqTrainMIHitFA] = gpfaAnalysis(Spikes.GPFA.MIHitFA.dat,6); %Run index
    close all
    [M1Dynamics(fileNum).neuralDynamics,waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);
    M1Dynamics(fileNum).filename = file(fileNum).name;
end
%% Trajectory Speed analysis
% Plot speed over time, highlighting different states
if ~exist('M1Dynamics','var') && ~exist('M2Dynamics','var')
    load('D:\TrajectoryDynamics\M1Dynamics.mat');
    load('D:\TrajectoryDynamics\M2Dynamics.mat');
end

dynamics = M2Dynamics;
dimension = 2;
speedTot = [];
for n = 1:length(dynamics)
speed_data = dynamics(n).neuralDynamics.hit.speed;
rt(n) = dynamics(n).neuralDynamics.mreactionTime  ;
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
speed_data = dynamics(n).neuralDynamics.miss.speed;
rt(n) = dynamics(n).neuralDynamics.mreactionTime  ;
speedTot{n} = squeeze(speed_data.speed(dimension,2:end,:));
end
speedTot = horzcat(speedTot{:});

plot(time(2:end),mean(speedTot,2),'color',colors(2,:),'linewidth',2),hold on
plot(time(2:end),mean(speedTot,2)+std(speedTot,[],2)/sqrt(size(speedTot,2)),'color',colors(2,:),'linewidth',2)
plot(time(2:end),mean(speedTot,2)-std(speedTot,[],2)/sqrt(size(speedTot,2)),'color',colors(2,:),'linewidth',2)

xline((75-75)*20, '--r', 'Cue');
xline((mean(rt)-80)*20, '--g', 'Movement Start');
xlabel('Time (s)');
ylabel('Average Speed');
box off,set(gca,'tickdir','out','fontsize',14),axis square,xlim([-500 1500]),ylim([0.01 0.05])

speedTot = [];
for n = 1:length(dynamics)
speed_data = dynamics(n).neuralDynamics.MIhit.speed;
rt(n) = dynamics(n).neuralDynamics.mreactionTime  ;
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
speed_data = dynamics(n).neuralDynamics.MIFA.speed;
rt(n) = dynamics(n).neuralDynamics.mreactionTime  ;
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
box off,set(gca,'tickdir','out','fontsize',14),axis square,xlim([-500 1500]),ylim([0.01 0.05])


%%
totalSpeedDimension = [];
for dimension = 1:6
speedHitCueMovement = [];
speedMissCueMovement = []; 
speedMIhitCueMovement = [];
speedFACueMovement = [];
for n = 1:length(dynamics)
speed_data = dynamics(n).neuralDynamics.hit.speed;
speedHitCueMovement{n} = squeeze(mean(speed_data.preMovement(dimension,:,:), [2 3]));

speed_data = dynamics(n).neuralDynamics.miss.speed;
speedMissCueMovement{n} = squeeze(mean(speed_data.preMovement(dimension,:,:), [2 3]));

speed_data = dynamics(n).neuralDynamics.MIhit.speed;
speedMIhitCueMovement{n} = squeeze(mean(speed_data.preMovement(dimension,:,:), [2 3]));

speed_data = dynamics(n).neuralDynamics.MIFA.speed;
speedFACueMovement{n} = squeeze(mean(speed_data.preMovement(dimension,:,:), [2 3]));

end
speedHitCueMovement = vertcat(speedHitCueMovement{:});
speedMissCueMovement = vertcat(speedMissCueMovement{:});
speedMIhitCueMovement = vertcat(speedMIhitCueMovement{:});
speedFACueMovement = vertcat(speedFACueMovement{:});


temp = nan(max([length(speedHitCueMovement),length(speedMissCueMovement),length(speedMIhitCueMovement),length(speedFACueMovement)]),4);
temp(1:length(speedHitCueMovement),1) = speedHitCueMovement;
temp(1:length(speedMissCueMovement),2) = speedMissCueMovement;
temp(1:length(speedMIhitCueMovement),3) = speedMIhitCueMovement;
temp(1:length(speedFACueMovement),4) = speedFACueMovement;

totalSpeedDimension{dimension} = temp;

end
%% Plot average speed for each state as a function of dimensions
colors = [0 0.4470 0.7410;0.75 0.75 0.75;0 0.4470 0.7410;190/255 30/255 45/255];

figure;
for n = 1:4
    dat = cellfun(@(x) horzcat(x(:,n)),totalSpeedDimension,'UniformOutput',false);
    dat = horzcat(dat{:});
    errorbar(1:6,nanmean(dat,1),nanstd(dat,[],1)./sqrt(size(dat,1)),'k.-','linewidth',2,'color',colors(n,:)),hold on
end
box off,set(gca,'tickdir','out','fontsize',14),axis square,xlim([0.5 6.5]),ylim([0 0.06])


[p, tbl, stats] = anova1(dat);

fprintf('ANOVA p-value: %f\n', p);

if p < 0.05
    c = multcompare(stats, 'Display', 'off');
    disp('Multiple Comparisons:');
    disp(c);
end


%%
figure;
customBarplot(temp,'scatter','off');
ylabel('Neural Trajectory Speed');
%title(['Average Trajectory Speed by Behavioral State for Dimension ' num2str(dimension)]);
box off,set(gca,'tickdir','out','fontsize',14),axis square,ylim([0 0.07])
%ylim([0 0.08])
% Perform statistical comparison (ANOVA)
[p, tbl, stats] = anova1(temp);

fprintf('ANOVA p-value: %f\n', p);

if p < 0.05
    c = multcompare(stats, 'Display', 'off');
    disp('Multiple Comparisons:');
    disp(c);
end
%%
neuralDynamics = neuralDynamics(1).neuralDynamics;
temp = nan(max([size(neuralDynamics.hit.s,2),size(neuralDynamics.miss.s,2),size(neuralDynamics.MIFA.s,2)]),3);
temp(1:size(neuralDynamics.hit.s,2),1) = neuralDynamics.hit.s(1,:);
temp(1:size(neuralDynamics.miss.s,2),2) = neuralDynamics.miss.s(1,:);
temp(1:size(neuralDynamics.MIFA.s,2),3) = neuralDynamics.MIFA.s(1,:);

figure,customBarplot(temp);
box off,set(gca,'tickdir','out','fontsize',14),axis square
ylabel('Trajectory Deviation')
[p,t,stats] = anova1(temp)
c = multcompare(stats)