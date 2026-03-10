plotFlag = 0;

for filenumber = 1:size(M2neuralDynamics,2)
    filenumber
    IntanBehaviour = M2neuralDynamics(filenumber).IntanBehaviour;
    neuralDynamics = M2neuralDynamics(filenumber).neuralDynamics;
    parameters = IntanBehaviour.parameters;
    rt = cell2mat(arrayfun(@(s) s.reactionTime, IntanBehaviour.cueHitTrace, 'UniformOutput', false));
    if mean(rt) < 0
        warning('Bad file')
        M2neuralDynamics(filenumber).skip = 1;
        continue
    else
        M2neuralDynamics(filenumber).skip = 0;
    end
    %% Matching Neural Dynamics to the IntanBehaviour file 
    % Removing disengage miss trials, hit trials with RT>1.5 and RT<0,
    % overlaping Hit and FA trials
    if ~isfield(IntanBehaviour,'cleanFlag')
        IntanBehaviour = cleanIntanBehaviour(IntanBehaviour);
    end
    if ~isfield(neuralDynamics,'matchFlag')
        % Removing the same trials as above from the neural trajectory data
        neuralDynamics = matchBehaviourNeuralDynamics(neuralDynamics,IntanBehaviour);
    end
    % Calculating neural trajectory speed
    neuralDynamics = calTrajSpeed(neuralDynamics,IntanBehaviour,IntanBehaviour.parameters);  
   
    %% Correlating behavior to the trajectory speed
    % Getting the peak speed from cue to MI in cue aligned hit traces 
    % Cue index will be cueIndex 
    % MI index will be MIIndex
    for i=1:size(IntanBehaviour.cueHitTrace,2)
        [neuralDynamics.hit.speed.cueMIPeakSpeed(i),neuralDynamics.hit.speed.cueMIPeakSpeedTime(i)] = max(neuralDynamics.hit.speed.speed(neuralDynamics.cueIndex:neuralDynamics.MIIndex(i),i));
        [neuralDynamics.hit.speed.cueEndPeakSpeed(i),neuralDynamics.hit.speed.cueEndPeakSpeedTime(i)] = max(neuralDynamics.hit.speed.speed(neuralDynamics.cueIndex:end,i));
        if neuralDynamics.RewardIndex(i) <= size(neuralDynamics.time,2)
            [neuralDynamics.hit.speed.cueRewardPeakSpeed(i),neuralDynamics.hit.speed.cueRewardPeakSpeedTime(i)] = max(neuralDynamics.hit.speed.speed(neuralDynamics.cueIndex:neuralDynamics.RewardIndex(i),i));
        else
            [neuralDynamics.hit.speed.cueRewardPeakSpeed(i),neuralDynamics.hit.speed.cueRewardPeakSpeedTime(i)] = max(neuralDynamics.hit.speed.speed(neuralDynamics.cueIndex:end,i));
        end
        [neuralDynamics.hitOnly.speed.cueMIPeakSpeed(i),neuralDynamics.hitOnly.speed.cueMIPeakSpeedTime(i)] = max(neuralDynamics.hitOnly.speed.speed(neuralDynamics.cueIndex:neuralDynamics.MIIndex(i),i));
        [neuralDynamics.hitOnly.speed.cueEndPeakSpeed(i),neuralDynamics.hitOnly.speed.cueEndPeakSpeedTime(i)] = max(neuralDynamics.hitOnly.speed.speed(neuralDynamics.cueIndex:end,i));
        if neuralDynamics.RewardIndex(i) <= size(neuralDynamics.time,2)
            [neuralDynamics.hitOnly.speed.cueRewardPeakSpeed(i),neuralDynamics.hitOnly.speed.cueRewardPeakSpeedTime(i)] = max(neuralDynamics.hitOnly.speed.speed(neuralDynamics.cueIndex:neuralDynamics.RewardIndex(i),i));
        else
            [neuralDynamics.hitOnly.speed.cueRewardPeakSpeed(i),neuralDynamics.hitOnly.speed.cueRewardPeakSpeedTime(i)] = max(neuralDynamics.hitOnly.speed.speed(neuralDynamics.cueIndex:end,i));
        end
    end
    
    % Getting the peak speed from cue to MI in MIHit traces 
    % Cue index will be cueIndex-(MIIndex-cueIndex) = 2*cueIndex-MIIndex
    % MI index will be cueIndex (midpoint)
    for i=1:size(IntanBehaviour.MIHitTrace,2)
        a = 2*neuralDynamics.cueIndex-neuralDynamics.MIIndex(i);
        if a < 1
            a = 1;
        end
        [neuralDynamics.MIhit.speed.cueMIPeakSpeed(i),neuralDynamics.MIhit.speed.cueMIPeakSpeedTime(i)] = max(neuralDynamics.hit.speed.speed(a:neuralDynamics.cueIndex,i));
    end
    
    % Movement time - Time from MI to the lever reaching the threshold
    neuralDynamics.movementTime = cell2mat(arrayfun(@(s) s.rewardIndex - s.LFPIndex(parameters.Fs*parameters.windowBeforePull), IntanBehaviour.MIHitTrace, 'UniformOutput', false));
    
    if plotFlag == 1
        % Fitting peak speed between cue and MI to reaction time 
        mdl = fitlm(neuralDynamics.hit.speed.cueMIPeakSpeed,IntanBehaviour.reactionTime)
        figure,plot(mdl);
        xlabel('Peak neural trajectory speed post cue');ylabel('Reaction Time');
        
        % Fitting peak speed between cue and MI to movement time 
        mdl = fitlm(neuralDynamics.MIhit.speed.cueMIPeakSpeed,neuralDynamics.movementTime)
        figure,plot(mdl);
        xlabel('Peak neural trajectory speed post cue');ylabel('Movement time');
        
        % Fitting time location of peak speed between cue and MI to reaction time 
        mdl = fitlm(neuralDynamics.hit.speed.cueMIPeakSpeedTime,IntanBehaviour.reactionTime)
        figure,plot(mdl);
        xlabel('Peak neural trajectory speed post cue');ylabel('Reaction Time');
    end
    
    %% Saving the IntanBehaviour and neuralDyanmics variable 
    M2neuralDynamics(filenumber).IntanBehaviour = IntanBehaviour;
    M2neuralDynamics(filenumber).neuralDynamics = neuralDynamics;
    %% Plotting
    if plotFlag == 1
        % Plotting Lever Trace and the neural trajectory speed - Cue alligned 
        figure;
        subplot(2,2,1);
        for i=1:size(neuralDynamics.hit.speed.speed,2)
            plot(neuralDynamics.time,neuralDynamics.hit.speed.speed(:,i),'Color',[0 0 0 0.1],'LineWidth',1.5);
            hold on;
        end
        plot(neuralDynamics.time,mean(neuralDynamics.hit.speed.speed,2),'Color',[1 0 0 1],'LineWidth',2);
        xline(1500,'--r','Cue','LabelVerticalAlignment','top');
        xline(1500+mean(IntanBehaviour.reactionTime,'all')*IntanBehaviour.parameters.Fs,'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
        xlabel('Time (in s)');box off;
        subplot(2,2,3);
        for i=1:size(IntanBehaviour.cueHitTrace,2)
            plot(IntanBehaviour.cueHitTrace(i).time,IntanBehaviour.cueHitTrace(i).trace,'Color',[0 0 0 0.1],'LineWidth',1.5);
            hold on;
        end
        plot(IntanBehaviour.cueHitTrace(1).time,mean(horzcat(IntanBehaviour.cueHitTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
        yline(IntanBehaviour.threshold,'--.b','Threshold','LabelHorizontalAlignment','left');
        xline(0,'--r','Cue','LabelVerticalAlignment','top');
        xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
        ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for Cue Hits');box off;
        
        subplot(2,2,2);
        for i=1:size(neuralDynamics.miss.speed.speed,2)
            plot(neuralDynamics.time,neuralDynamics.miss.speed.speed(:,i),'Color',[0 0 0 0.1],'LineWidth',1.5);
            hold on;
        end
        plot(neuralDynamics.time,mean(neuralDynamics.miss.speed.speed,2),'Color',[1 0 0 1],'LineWidth',2);
        xline(1500,'--r','Cue','LabelVerticalAlignment','top');
        xline(1500+mean(IntanBehaviour.reactionTime,'all')*IntanBehaviour.parameters.Fs,'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
        xlabel('Time (in s)');box off;
        subplot(2,2,4);
        for i=1:size(IntanBehaviour.cueMissTrace,2)
            plot(IntanBehaviour.cueMissTrace(i).time,IntanBehaviour.cueMissTrace(i).trace,'Color',[0 0 0 0.1],'LineWidth',1.5);
            hold on;
        end
        plot(IntanBehaviour.cueMissTrace(1).time,mean(horzcat(IntanBehaviour.cueMissTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
        yline(IntanBehaviour.threshold,'--.b','Threshold','LabelHorizontalAlignment','left');
        xline(0,'--r','Cue','LabelVerticalAlignment','top');
        xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
        ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for Cue Misss');box off;
        
        % Plotting Lever Trace and the neural trajectory speed - MI alligned 
        figure;
        subplot(2,2,1);
        for i=1:size(neuralDynamics.MIhit.speed.speed,2)
            plot(neuralDynamics.time,neuralDynamics.MIhit.speed.speed(:,i),'Color',[0 0 0 0.1],'LineWidth',1.5);
            hold on;
        end
        plot(neuralDynamics.time,mean(neuralDynamics.MIhit.speed.speed,2),'Color',[1 0 0 1],'LineWidth',2);
        xline(1500,'--r','MI','LabelVerticalAlignment','top');
        xlabel('Time (in s)');box off;ylabel('Trajectory speed');
        subplot(2,2,3);
        for i=1:size(IntanBehaviour.MIHitTrace,2)
            plot(IntanBehaviour.MIHitTrace(i).time,IntanBehaviour.MIHitTrace(i).trace,'Color',[0 0 0 0.1],'LineWidth',1.5);
            hold on;
        end
        plot(IntanBehaviour.MIHitTrace(1).time,mean(horzcat(IntanBehaviour.MIHitTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
        yline(IntanBehaviour.MIcutoffHit,'--.b','MI Threshold','LabelHorizontalAlignment','left'); 
        yline(IntanBehaviour.threshold,'--.b','Reward Threshold','LabelHorizontalAlignment','left'); 
        xline(0,'--r','MI','LabelVerticalAlignment','top');
        ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for MI Hits');box off;
        
        subplot(2,2,2);
        for i=1:size(neuralDynamics.MIFA.speed.speed,2)
            plot(neuralDynamics.time,neuralDynamics.MIFA.speed.speed(:,i),'Color',[0 0 0 0.1],'LineWidth',1.5);
            hold on;
        end
        plot(neuralDynamics.time,mean(neuralDynamics.MIFA.speed.speed,2),'Color',[1 0 0 1],'LineWidth',2);
        xline(1500,'--r','MI','LabelVerticalAlignment','top');
        xlabel('Time (in s)');box off;ylabel('Trajectory speed');
        subplot(2,2,4);
        for i=1:size(IntanBehaviour.MIFATrace,2)
            plot(IntanBehaviour.MIFATrace(i).time,IntanBehaviour.MIFATrace(i).trace,'Color',[0 0 0 0.1],'LineWidth',1.5);
            hold on;
        end
        plot(IntanBehaviour.MIFATrace(1).time,mean(horzcat(IntanBehaviour.MIFATrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
        yline(IntanBehaviour.MIcutoffHit,'--.b','MI Threshold','LabelHorizontalAlignment','left'); 
        yline(IntanBehaviour.threshold,'--.b','Reward Threshold','LabelHorizontalAlignment','left'); 
        xline(0,'--r','MI','LabelVerticalAlignment','top');
        ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for MI FA');box off;
    end
end

%% Pooling data from all animals and sessions 

% SHOULD I NORMALIZE FOR EACH MOUSE?
PooledData.hit.speed.cueMIPeakSpeedTime = [];
PooledData.hit.speed.cueMIPeakSpeed = [];
PooledData.hit.speed.cueEndPeakSpeedTime = [];
PooledData.hit.speed.cueEndPeakSpeed = [];
PooledData.hit.speed.cueRewardPeakSpeedTime = [];
PooledData.hit.speed.cueRewardPeakSpeed = [];
PooledData.RT = [];
for i = 1:size(M2neuralDynamics,2)
    try
    if M2neuralDynamics(i).skip == 1
        continue
    else
        PooledData.RT = [PooledData.RT ,  M2neuralDynamics(i).IntanBehaviour.reactionTime];
        PooledData.hit.speed.cueMIPeakSpeedTime = ([PooledData.hit.speed.cueMIPeakSpeedTime , M2neuralDynamics(i).neuralDynamics.hit.speed.cueMIPeakSpeedTime ]);
        PooledData.hit.speed.cueMIPeakSpeed = ([PooledData.hit.speed.cueMIPeakSpeed , M2neuralDynamics(i).neuralDynamics.hit.speed.cueMIPeakSpeed ]);
        PooledData.hit.speed.cueEndPeakSpeedTime = ([PooledData.hit.speed.cueEndPeakSpeedTime , M2neuralDynamics(i).neuralDynamics.hit.speed.cueEndPeakSpeedTime ]);
        PooledData.hit.speed.cueEndPeakSpeed = ([PooledData.hit.speed.cueEndPeakSpeed , M2neuralDynamics(i).neuralDynamics.hit.speed.cueEndPeakSpeed ]);
        PooledData.hit.speed.cueRewardPeakSpeedTime = ([PooledData.hit.speed.cueRewardPeakSpeedTime , M2neuralDynamics(i).neuralDynamics.hit.speed.cueRewardPeakSpeedTime ]);
        PooledData.hit.speed.cueRewardPeakSpeed = ([PooledData.hit.speed.cueRewardPeakSpeed , M2neuralDynamics(i).neuralDynamics.hit.speed.cueRewardPeakSpeed ]);
    end
    catch
        continue
    end
end
% Fitting time location of peak speed between cue and MI to reaction time 
neuralTrajectoryTs = 20; % in ms
mdl = fitlm(neuralTrajectoryTs*PooledData.hit.speed.cueMIPeakSpeedTime,PooledData.RT)
figure,plot(mdl);
xlabel('Time of peak neural trajectory speed from cue to MI');ylabel('Reaction Time');

% Fitting peak speed between cue and MI to reaction time 
mdl = fitlm(PooledData.hit.speed.cueMIPeakSpeed,PooledData.RT)
figure,plot(mdl);
xlabel('Peak neural trajectory speed from cue');ylabel('Reaction Time');

% Removing outlier in RT 
% Fitting time location of peak speed between cue and MI to reaction time 
[ ~,outlierIndex] = rmoutliers(PooledData.RT);
PooledData.hit.speed.cueMIPeakSpeedRMOutlier = neuralTrajectoryTs*PooledData.hit.speed.cueMIPeakSpeed;
PooledData.hit.speed.cueMIPeakSpeedRMOutlier(outlierIndex) = [];
PooledData.RTRMOutlier = PooledData.RT; PooledData.RTRMOutlier(outlierIndex)=[];
mdl = fitlm(PooledData.hit.speed.cueMIPeakSpeedRMOutlier,PooledData.RTRMOutlier)
figure,plot(mdl);
xlabel('Peak neural trajectory speed from cue');ylabel('Reaction Time');

% Fitting time location of peak speed between cue and end to reaction time 
neuralTrajectoryTs = 20; % in ms
mdl = fitlm(neuralTrajectoryTs*PooledData.hit.speed.cueEndPeakSpeedTime,PooledData.RT)
figure,plot(mdl);
xlabel('Time of peak neural trajectory speed from cue to End');ylabel('Reaction Time');

% Fitting time location of peak speed between cue and reward to reaction time 
neuralTrajectoryTs = 20; % in ms
mdl = fitlm(neuralTrajectoryTs*PooledData.hit.speed.cueRewardPeakSpeedTime,PooledData.RT)
figure,plot(mdl);
xlabel('Time of peak neural trajectory speed from cue to Reward');ylabel('Reaction Time');

%% z-Scoring plotting 
PooledData.time = M2neuralDynamics(1).neuralDynamics.time;
PooledData.hit.speed.zspeed = [];
PooledData.miss.speed.zspeed = [];
PooledData.hitOnly.speed.zspeed = [];
PooledData.MIhit.speed.zspeed = [];
PooledData.MIFA.speed.zspeed = [];

for i = 1:size(M2neuralDynamics,2)
    if M2neuralDynamics(i).skip == 1
        continue
    else
        PooledData.hit.speed.zspeed = ([PooledData.hit.speed.zspeed  M2neuralDynamics(i).neuralDynamics.hit.speed.zspeed]);
        PooledData.hitOnly.speed.zspeed = ([PooledData.hitOnly.speed.zspeed  M2neuralDynamics(i).neuralDynamics.hitOnly.speed.zspeed]);
        PooledData.miss.speed.zspeed = ([PooledData.miss.speed.zspeed  M2neuralDynamics(i).neuralDynamics.miss.speed.zspeed]);
        PooledData.MIhit.speed.zspeed = ([PooledData.MIhit.speed.zspeed  M2neuralDynamics(i).neuralDynamics.MIhit.speed.zspeed]);
        PooledData.MIFA.speed.zspeed = ([PooledData.MIFA.speed.zspeed  M2neuralDynamics(i).neuralDynamics.MIFA.speed.zspeed]);
    end
end

% Plotting z-scored trajectory speeds 
figure();
subplot(3,1,1);
h1=plot(PooledData.time,squeeze(mean(PooledData.hit.speed.zspeed,2)),'Color', [0 0.1 0.8],'LineWidth',2); hold on;
plot(PooledData.time,squeeze(mean(PooledData.hit.speed.zspeed,2)) - (squeeze(std(PooledData.hit.speed.zspeed,0,2)))/sqrt(size(PooledData.hit.speed.zspeed,2)) ,'Color', [0 0.1 0.8 0.4],'LineWidth',2);
plot(PooledData.time,squeeze(mean(PooledData.hit.speed.zspeed,2)) + (squeeze(std(PooledData.hit.speed.zspeed,0,2)))/sqrt(size(PooledData.hit.speed.zspeed,2)) ,'Color', [0 0.1 0.8 0.4],'LineWidth',2);
h2=plot(PooledData.time,squeeze(mean(PooledData.miss.speed.zspeed,2)),'Color', [0.5 0.5 0.5],'LineWidth',2); hold on;
plot(PooledData.time,squeeze(mean(PooledData.miss.speed.zspeed,2)) - (squeeze(std(PooledData.miss.speed.zspeed,0,2)))/sqrt(size(PooledData.miss.speed.zspeed,2)) ,'Color', [0.5 0.5 0.5 0.4],'LineWidth',2);
plot(PooledData.time,squeeze(mean(PooledData.miss.speed.zspeed,2)) + (squeeze(std(PooledData.miss.speed.zspeed,0,2)))/sqrt(size(PooledData.miss.speed.zspeed,2)) ,'Color', [0.5 0.5 0.5 0.4],'LineWidth',2);
xline(1501,'--r','Cue');%xline(1500+parameters.Fs*mean(IntanBehaviour.reactionTime,'all'),'--r','RT');
xlabel('Time (ms)'); ylabel('Neural Trajectory Speed (z-scored)');
legend([h1 h2],'Hit','Miss','Location','best'); %ylim([5 15]);
xlim([0 3000]);box off;set(gca,'TickDir','out','fontsize',14');
drawnow;

subplot(3,1,2);
h1=plot(PooledData.time,squeeze(mean(PooledData.MIhit.speed.zspeed,2)),'Color', [0 0.1 0.8],'LineWidth',2); hold on;
plot(PooledData.time,squeeze(mean(PooledData.MIhit.speed.zspeed,2)) - (squeeze(std(PooledData.MIhit.speed.zspeed,0,2)))/sqrt(size(PooledData.MIhit.speed.zspeed,2)) ,'Color', [0 0.1 0.8 0.4],'LineWidth',2);
plot(PooledData.time,squeeze(mean(PooledData.MIhit.speed.zspeed,2)) + (squeeze(std(PooledData.MIhit.speed.zspeed,0,2)))/sqrt(size(PooledData.MIhit.speed.zspeed,2)) ,'Color', [0 0.1 0.8 0.4],'LineWidth',2);
h2=plot(PooledData.time,squeeze(mean(PooledData.MIFA.speed.zspeed,2)),'Color', [0.9 0.1 0.1],'LineWidth',2); hold on;
plot(PooledData.time,squeeze(mean(PooledData.MIFA.speed.zspeed,2)) - (squeeze(std(PooledData.MIFA.speed.zspeed,0,2)))/sqrt(size(PooledData.MIFA.speed.zspeed,2)) ,'Color', [0.9 0.1 0.1 0.4],'LineWidth',2);
plot(PooledData.time,squeeze(mean(PooledData.MIFA.speed.zspeed,2)) + (squeeze(std(PooledData.MIFA.speed.zspeed,0,2)))/sqrt(size(PooledData.MIFA.speed.zspeed,2)) ,'Color', [0.9 0.1 0.1 0.4],'LineWidth',2);
xline(1501,'--r','MI');%xline(1500+parameters.Fs*mean(IntanBehaviour.reactionTime,'all'),'--r','RT');
xlabel('Time (ms)'); ylabel('Neural Trajectory Speed (z-scored)');
legend([h1 h2],'MI Hit','MI FA','Location','best'); %ylim([5 15]);
xlim([0 3000]);box off;set(gca,'TickDir','out','fontsize',14');
drawnow;

subplot(3,1,3);
h1=plot(PooledData.time,squeeze(mean(PooledData.hitOnly.speed.zspeed,2)),'Color', [0 0.1 0.8],'LineWidth',2); hold on;
plot(PooledData.time,squeeze(mean(PooledData.hitOnly.speed.zspeed,2)) - (squeeze(std(PooledData.hitOnly.speed.zspeed,0,2)))/sqrt(size(PooledData.hitOnly.speed.zspeed,2)) ,'Color', [0 0.1 0.8 0.4],'LineWidth',2);
plot(PooledData.time,squeeze(mean(PooledData.hitOnly.speed.zspeed,2)) + (squeeze(std(PooledData.hitOnly.speed.zspeed,0,2)))/sqrt(size(PooledData.hitOnly.speed.zspeed,2)) ,'Color', [0 0.1 0.8 0.4],'LineWidth',2);
xline(1501,'--r','Cue');%xline(1500+parameters.Fs*mean(IntanBehaviour.reactionTime,'all'),'--r','RT');
xlabel('Time (ms)'); ylabel('Neural Trajectory Speed (z-scored)');
legend([h1],'Hit Only','Location','best'); %ylim([5 15]);
xlim([0 3000]);box off;set(gca,'TickDir','out','fontsize',14');
drawnow;

%% LOCAL FUNCTIONs
function [neuralDynamics] = matchBehaviourNeuralDynamics(neuralDynamics,IntanBehaviour)
% Rejecting neural trajectories according to IntanBehaviour 
% Hit
neuralDynamics.hit.X(:,:,IntanBehaviour.rejectTrials.Hit) = [];
% Miss
neuralDynamics.miss.X(:,:,IntanBehaviour.rejectTrials.cueMiss) = [];
% MIHit
neuralDynamics.MIhit.X(:,:,IntanBehaviour.rejectTrials.Hit) = [];
% MIFA
neuralDynamics.MIFA.X(:,:,IntanBehaviour.rejectTrials.FA) = [];

if isfield(neuralDynamics,'hitOnly')
    % Hit only trajectory
    neuralDynamics.hitOnly.X(:,:,IntanBehaviour.rejectTrials.Hit) = [];
end

neuralDynamics.matchFlag = 1;
end

function [neuralDynamics] = calTrajSpeed(neuralDynamics,IntanBehaviour,parameters)

% Calculating neural trajectory speed 
neuralDynamics.nPCA = size(neuralDynamics.hit.X,1);
ntime = ((parameters.windowBeforeCue+parameters.windowAfterCue)*parameters.Fs)/size(neuralDynamics.hit.X,2);  
neuralDynamics.time = ntime:ntime:((parameters.windowBeforeCue+parameters.windowAfterCue)*parameters.Fs);

neuralDynamics.preCueIndex = max(find(neuralDynamics.time<parameters.windowBeforeCue*parameters.Fs));
neuralDynamics.cueIndex = neuralDynamics.preCueIndex+1;
RT = (parameters.windowBeforeCue + cell2mat(arrayfun(@(s) s.reactionTime, IntanBehaviour.cueHitTrace, 'UniformOutput', false)))*parameters.Fs;
neuralDynamics.MIIndex = round(RT./ntime);
RewardT = (parameters.windowBeforeCue*parameters.Fs + cell2mat(arrayfun(@(s) (s.rewardIndex - s.LFPIndex(ceil(end/2))), IntanBehaviour.cueHitTrace, 'UniformOutput', false)));
neuralDynamics.RewardIndex = round(RewardT./ntime);
neuralDynamics.cueIndexMI = neuralDynamics.cueIndex - (floor(mean(neuralDynamics.MIIndex))-neuralDynamics.cueIndex);
% neuralDynamics.timeSpeed = neuralDynamics.time(1:end-1);
%% For hit 
neuralDynamics.hit.speed.vel = diff(neuralDynamics.hit.X,1,2)/ntime;
neuralDynamics.hit.speed.speed = squeeze(sqrt(sum(neuralDynamics.hit.speed.vel.^2,1)));
% pad zeros in front 
neuralDynamics.hit.speed.speed = [zeros(1,size(neuralDynamics.hit.speed.speed,2));neuralDynamics.hit.speed.speed];
% neuralDynamics.hit.speed.preCueSpeed = neuralDynamics.hit.speed.speed(1:preCueIndex,:);

%% For miss
neuralDynamics.miss.speed.vel = diff(neuralDynamics.miss.X,1,2)/ntime;
neuralDynamics.miss.speed.speed = squeeze(sqrt(sum(neuralDynamics.miss.speed.vel.^2,1)));
% pad zeros in front 
neuralDynamics.miss.speed.speed = [zeros(1,size(neuralDynamics.miss.speed.speed,2));neuralDynamics.miss.speed.speed];
% neuralDynamics.miss.speed.preCueSpeed = neuralDynamics.miss.speed.speed(1:preCueIndex,:);

%% ZScoring Hits vs Misses
muhitmiss = mean([neuralDynamics.hit.speed.speed(2:neuralDynamics.preCueIndex,:) neuralDynamics.miss.speed.speed(2:neuralDynamics.preCueIndex,:)],'all');
stdhitmiss = std([neuralDynamics.hit.speed.speed(2:neuralDynamics.preCueIndex,:) neuralDynamics.miss.speed.speed(2:neuralDynamics.preCueIndex,:)],0,'all');
neuralDynamics.hit.speed.zspeed = (neuralDynamics.hit.speed.speed - muhitmiss)/stdhitmiss;
neuralDynamics.miss.speed.zspeed = (neuralDynamics.miss.speed.speed - muhitmiss)/stdhitmiss;
neuralDynamics.hit.speed.zspeed(1,:) = 0;
neuralDynamics.miss.speed.zspeed(1,:) = 0;
%% For MIhit 
neuralDynamics.MIhit.speed.vel = diff(neuralDynamics.MIhit.X,1,2)/ntime;
neuralDynamics.MIhit.speed.speed = squeeze(sqrt(sum(neuralDynamics.MIhit.speed.vel.^2,1)));
% pad zeros in front 
neuralDynamics.MIhit.speed.speed = [zeros(1,size(neuralDynamics.MIhit.speed.speed,2));neuralDynamics.MIhit.speed.speed];

%% For MIFA
neuralDynamics.MIFA.speed.vel = diff(neuralDynamics.MIFA.X,1,2)/ntime;
neuralDynamics.MIFA.speed.speed = squeeze(sqrt(sum(neuralDynamics.MIFA.speed.vel.^2,1)));
% pad zeros in front 
neuralDynamics.MIFA.speed.speed = [zeros(1,size(neuralDynamics.MIFA.speed.speed,2));neuralDynamics.MIFA.speed.speed];

%% ZScoring MIHits vs MIFAs
muhitfa = mean([neuralDynamics.MIhit.speed.speed(2:neuralDynamics.cueIndexMI,:) neuralDynamics.MIFA.speed.speed(2:neuralDynamics.cueIndexMI,:)],'all');
stdhitfa = std([neuralDynamics.MIhit.speed.speed(2:neuralDynamics.cueIndexMI,:) neuralDynamics.MIFA.speed.speed(2:neuralDynamics.cueIndexMI,:)],0,'all');
neuralDynamics.MIhit.speed.zspeed = (neuralDynamics.MIhit.speed.speed - muhitfa)/stdhitfa;
neuralDynamics.MIFA.speed.zspeed = (neuralDynamics.MIFA.speed.speed - muhitfa)/stdhitfa;
neuralDynamics.MIhit.speed.zspeed(1,:) = 0;
neuralDynamics.MIFA.speed.zspeed(1,:) = 0;

%% for hit only trajectory
if isfield(neuralDynamics,'hitOnly')
    neuralDynamics.hitOnly.speed.vel = diff(neuralDynamics.hitOnly.X,1,2)/ntime;
    neuralDynamics.hitOnly.speed.speed = squeeze(sqrt(sum(neuralDynamics.hitOnly.speed.vel.^2,1)));
    muhitonly = mean(neuralDynamics.hitOnly.speed.speed(2:neuralDynamics.preCueIndex,:),'all');
    stdhitonly = std(neuralDynamics.hitOnly.speed.speed(2:neuralDynamics.preCueIndex,:),0,'all');
    neuralDynamics.hitOnly.speed.zspeed = (neuralDynamics.hitOnly.speed.speed-muhitonly)/stdhitonly;
    % pad zeros in front 
    neuralDynamics.hitOnly.speed.speed = [zeros(1,size(neuralDynamics.hitOnly.speed.speed,2));neuralDynamics.hitOnly.speed.speed];
    neuralDynamics.hitOnly.speed.zspeed = [zeros(1,size(neuralDynamics.hitOnly.speed.zspeed,2));neuralDynamics.hitOnly.speed.zspeed];
end


% neuralDynamics.hit.speed.preCueSpeed = neuralDynamics.hit.speed.speed(1:preCueIndex,:);
% figure;
% for i=1:size(neuralDynamics.hit.speed,2)
%     plot(neuralDynamics.timeSpeed,neuralDynamics.hit.speed(:,i),'Color',[0 0 0 0.1],'LineWidth',1.5);
%     hold on;
% end
% plot(neuralDynamics.timeSpeed,mean(neuralDynamics.hit.speed,2),'Color',[1 0 0 1],'LineWidth',2);


end