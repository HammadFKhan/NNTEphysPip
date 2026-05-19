clear
%fpath = 'F:\LeverTask\Ephys\Analysis\M2Spikes';
%fpath = 'F:\LeverTask\Ephys\Analysis\spksPooledwFA'
%fpath = 'D:\M1_GSP';
% fpath = 'D:\M2SpikeData';
% file = dir(fullfile(fpath,'*.mat'));
clear
redo = 1;
if redo==1
    files = dir(fullfile('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\M1_GSP','*.mat'));
    M2neuralDynamics = struct();
    for fileNum = 1:length(files)
        disp(['File number: ' num2str(fileNum)])
        load(fullfile(files(fileNum).folder,files(fileNum).name))
        try
        Spikes = makeSpikeGPFA(Spikes);
        Spikes.GPFA.HitMiss.dat = [Spikes.GPFA.hit.dat,Spikes.GPFA.miss.dat];
        for n = length(Spikes.GPFA.hit.dat)+1:length(Spikes.GPFA.HitMiss.dat) %fix trials
            Spikes.GPFA.HitMiss.dat(n).trialId = n;
        end
        Spikes.GPFA.MIHitFA.dat = [Spikes.GPFA.MIHit.dat,Spikes.GPFA.MIFA.dat];
        for n = length(IntanBehaviour.MIHitTrace)+1:length(Spikes.GPFA.MIHitFA.dat) %fix trials
            Spikes.GPFA.MIHitFA.dat(n).trialId = n;
        end
        addpath(genpath('C:\Users\khan332\Documents\GitHub\NeuralTraj'));
        addpath(genpath('mat_results'));
        if exist('mat_results','dir'),rmdir('mat_results','s'),end
        [Spikes.GPFA.resultHit,Spikes.GPFA.seqTrainHit] = gpfaAnalysis(Spikes.GPFA.hit.dat,1); %Run index
        [Spikes.GPFA.resultMiss,Spikes.GPFA.seqTrainMiss] = gpfaAnalysis(Spikes.GPFA.miss.dat,2); %Run index
        [Spikes.GPFA.resultMIHit,Spikes.GPFA.seqTrainMIHit] = gpfaAnalysis(Spikes.GPFA.MIHit.dat,3); %Run index
        [Spikes.GPFA.resultMIFA,Spikes.GPFA.seqTrainMIFA] = gpfaAnalysis(Spikes.GPFA.MIFA.dat,4); %Run index
        [Spikes.GPFA.resultHitMiss,Spikes.GPFA.seqTrainHitMiss] = gpfaAnalysis(Spikes.GPFA.HitMiss.dat,5); %Run index
        [Spikes.GPFA.resultMIHitFA,Spikes.GPFA.seqTrainMIHitFA] = gpfaAnalysis(Spikes.GPFA.MIHitFA.dat,6); %Run index
        M2neuralDynamics(fileNum).fname = files(fileNum).name;
        M2neuralDynamics(fileNum).IntanBehaviour = IntanBehaviour;
        [M2neuralDynamics(fileNum).neuralDynamics,M1waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);
        catch ME
            disp('Error getting neural dynamics')
            continue
        end
        close all
    end
    fpath = 'Y:\Hammad\Ephys\LeverTask\Data_for_Figures\TrajectoryDynamics';
    sessionName = [fpath,'\','M2DynamicsPooledRTv2.mat'];
    save(sessionName,"M2neuralDynamics","fileNum","files","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
    clear
    disp('Loading processed data...')
    load('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\TrajectoryDynamics\M2DynamicsPooledRTv2.mat')
else
    fprintf('Loading processed data...')
    load('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\TrajectoryDynamics\M2DynamicsPooledRTv2.mat')
    fprintf('done\n')
end
%% Trajectory Speed analysis
% Plot speed over time, highlighting different states
if ~exist('M1Dynamics','var') && ~exist('M2Dynamics','var')
    load('D:\TrajectoryDynamics\M1Dynamics.mat');
    load('D:\TrajectoryDynamics\M2Dynamics.mat');
end
%%
dynamics = M2neuralDynamics;
dimension = 1;
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
box off,set(gca,'tickdir','out','fontsize',14),axis square,xlim([-500 1500]),ylim([0.015 0.05])


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
%% Neural Trajectory difference
dynamics = M1neuralDynamics;
dimension = 1;
simTot = [];
for n = 1:length(dynamics)
sim_data = dynamics(n).neuralDynamics_unbalanced.neuralDiffhitmiss;  
simTot{n} = squeeze(sim_data(:,dimension));
end
simTot = horzcat(simTot{:});


colors = [0 0.4470 0.7410;0.75 0.75 0.75;190/255 30/255 45/255];
time = -1499:20:1500;
figure;
plot(time(2:end),mean(simTot,2),'color',colors(1,:),'linewidth',2),hold on
plot(time(2:end),mean(simTot,2)+std(simTot,[],2)/sqrt(size(simTot,2)),'color',colors(1,:),'linewidth',2)
plot(time(2:end),mean(simTot,2)-std(simTot,[],2)/sqrt(size(simTot,2)),'color',colors(1,:),'linewidth',2)
hold on;

dimension = 1;
simTot = [];
for n = 1:length(dynamics)
sim_data = dynamics(n).neuralDynamics_rebalanced.neuralDiffhitmiss;  
simTot{n} = squeeze(sim_data(:,dimension));
end
simTot = horzcat(simTot{:});


colors = [0 0.4470 0.7410;0.75 0.75 0.75;190/255 30/255 45/255];
time = -1499:20:1500;
plot(time(2:end),mean(simTot,2),'color',colors(2,:),'linewidth',2),hold on
plot(time(2:end),mean(simTot,2)+std(simTot,[],2)/sqrt(size(simTot,2)),'color',colors(2,:),'linewidth',2)
plot(time(2:end),mean(simTot,2)-std(simTot,[],2)/sqrt(size(simTot,2)),'color',colors(2,:),'linewidth',2)
hold on;

xlabel('Time (s)');
ylabel('Average Speed');
box off,set(gca,'tickdir','out','fontsize',14),axis square,xlim([-500 1500])%,ylim([-0.01 0.03
%%
dynamics = M1neuralDynamics;
dimension = 1;
simTot = [];
for n = 1:length(dynamics)
sim_data = dynamics(n).neuralDynamics_unbalanced.neuralDiffMI;  
simTot{n} = squeeze(sim_data(:,dimension));
end
simTot = horzcat(simTot{:});


colors = [0 0.4470 0.7410;0.75 0.75 0.75;190/255 30/255 45/255];
time = -1499:20:1500;
figure;
plot(time(2:end),mean(simTot,2),'color',colors(1,:),'linewidth',2),hold on
plot(time(2:end),mean(simTot,2)+std(simTot,[],2)/sqrt(size(simTot,2)),'color',colors(1,:),'linewidth',2)
plot(time(2:end),mean(simTot,2)-std(simTot,[],2)/sqrt(size(simTot,2)),'color',colors(1,:),'linewidth',2)
hold on;

simTot = [];
for n = 1:length(dynamics)
sim_data = dynamics(n).neuralDynamics_rebalanced.neuralDiffMI;  
simTot{n} = squeeze(sim_data(:,dimension));
end
simTot = horzcat(simTot{:});


colors = [0 0.4470 0.7410;0.75 0.75 0.75;190/255 30/255 45/255];
time = -1499:20:1500;
plot(time(2:end),mean(simTot,2),'color',colors(2,:),'linewidth',2),hold on
plot(time(2:end),mean(simTot,2)+std(simTot,[],2)/sqrt(size(simTot,2)),'color',colors(2,:),'linewidth',2)
plot(time(2:end),mean(simTot,2)-std(simTot,[],2)/sqrt(size(simTot,2)),'color',colors(2,:),'linewidth',2)
hold on;

xlabel('Time (s)');
ylabel('Average Speed');
box off,set(gca,'tickdir','out','fontsize',14),axis square,xlim([-500 1500])%,ylim([-0.01 0.03])
