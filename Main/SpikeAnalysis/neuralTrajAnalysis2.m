function [neuralDynamics,waveDynamics] = neuralTrajAnalysis2(Spikes,Waves,Behaviour)
%% Take orthoganal latent dimensions and do statistics across trials
X = arrayfun(@(x) vertcat(x.xorth),Spikes.GPFA.seqTrainHit,'UniformOutput',false);
X = horzcat(X{:});
neuralTrajHit = reshape(X,size(X,1),Spikes.GPFA.seqTrainHit(1).T,[]);

X = arrayfun(@(x) vertcat(x.xorth),Spikes.GPFA.seqTrainMiss,'UniformOutput',false);
X = horzcat(X{:});
neuralTrajMiss = reshape(X,size(X,1),Spikes.GPFA.seqTrainMiss(1).T,[]);

X = arrayfun(@(x) vertcat(x.xorth),Spikes.GPFA.seqTrainMIHit,'UniformOutput',false);
X = horzcat(X{:});
neuralTrajMIHit = reshape(X,size(X,1),Spikes.GPFA.seqTrainMIHit(1).T,[]);

X = arrayfun(@(x) vertcat(x.xorth),Spikes.GPFA.seqTrainMIFA,'UniformOutput',false);
X = horzcat(X{:});
neuralTrajMIFA = reshape(X,size(X,1),Spikes.GPFA.seqTrainMIFA(1).T,[]);

X = arrayfun(@(x) vertcat(x.xorth),Spikes.GPFA.seqTrainHitMiss,'UniformOutput',false);
X = horzcat(X{:});
neuralTrajHitMiss = reshape(X,size(X,1),Spikes.GPFA.seqTrainHitMiss(1).T,[]);

X = arrayfun(@(x) vertcat(x.xorth),Spikes.GPFA.seqTrainMIHitFA,'UniformOutput',false);
X = horzcat(X{:});
neuralTrajMIHitFA = reshape(X,size(X,1),Spikes.GPFA.seqTrainMIHitFA(1).T,[]);

% X = arrayfun(@(x) vertcat(x.xorth),Spikes.GPFA.seqTrainBaselineOpto,'UniformOutput',false);
% X = horzcat(X{:});
% neuralTrajBaselineOpto = reshape(X,size(X,1),Spikes.GPFA.seqTrainBaselineOpto(1).T,[]);

%% Calculate average trajectories and divergence based on trial difference
% local function call for meaning based on combined PCA of trial conditions
dimNum = 15; %nu,ber of dimensions to take
X = neuralTrajHitMiss;
hittrials = 1:length(Behaviour.cueHitTrace);
misstrials = length(Behaviour.cueHitTrace)+1:size(X,3);

[neuralDynamics.hit.r,neuralDynamics.hit.s,neuralDynamics.hit.stab,neuralDynamics.hit.X] = getMeanTraj(X,hittrials,dimNum); %trajectory variable and predefined conditional trial indexes
neuralDynamics.hit.speed = speedTraj(X,hittrials,6,Behaviour);

[neuralDynamics.miss.r,neuralDynamics.miss.s,neuralDynamics.miss.stab,neuralDynamics.miss.X] = getMeanTraj(X,misstrials,dimNum); %trajectory variable and predefined conditional trial indexes
neuralDynamics.miss.speed = speedTraj(X,misstrials,6,Behaviour);

% Calculate differences in trajectories r'c(t)/||r'c(t)||
[neuralDynamics.neuralSimhitmiss,neuralDynamics.neuralDiffhitmiss,neuralDynamics.hit.rprime,neuralDynamics.miss.rprime] = neuralTrajDiff(neuralDynamics.hit.r',neuralDynamics.miss.r');

%% Do just for hits
X = neuralTrajHit;
hittrials = 1:length(Behaviour.cueHitTrace);
[neuralDynamics.hit.r,neuralDynamics.hitOnly.s,neuralDynamics.hitOnly.stab,neuralDynamics.hitOnly.X] = getMeanTraj(X,hittrials,dimNum); %trajectory variable and predefined conditional trial indexes
neuralDynamics.hitOnly.speed = speedTraj(X,hittrials,6,Behaviour);
%% Now do analysis for MI hit vs FA
X = neuralTrajMIHitFA;
MIHittrials = 1:length(Behaviour.MIHitTrace);
MIFAtrials = length(Behaviour.MIHitTrace)+1:size(X,3);

[neuralDynamics.MIhit.r,neuralDynamics.MIhit.s,neuralDynamics.MIhit.stab,neuralDynamics.MIhit.X] = getMeanTraj(X,MIHittrials,dimNum); %trajectory variable and predefined conditional trial indexes
neuralDynamics.MIhit.speed = speedTraj(X,MIHittrials,6,Behaviour);


[neuralDynamics.MIFA.r,neuralDynamics.MIFA.s,neuralDynamics.MIFA.stab,neuralDynamics.MIFA.X] = getMeanTraj(X,MIFAtrials,dimNum); %trajectory variable and predefined conditional trial indexes
neuralDynamics.MIFA.speed = speedTraj(X,MIFAtrials,6,Behaviour);

% Calculate differences in trajectories r'c(t)/||r'c(t)||
[neuralDynamics.neuralSimMI,neuralDynamics.neuralDiffMI,neuralDynamics.MIhit.rprime,neuralDynamics.MIFA.rprime] = neuralTrajDiff(neuralDynamics.MIhit.r',neuralDynamics.MIFA.r');
%%
% % rhminit = abs(rh(1,:)-rm(1,:));
% % rm = rm-rhminit;
t = linspace(-Behaviour.parameters.windowBeforeCue,Behaviour.parameters.windowAfterCue,size(neuralDynamics.hit.r,2));
% Find important indices in array
neuralDynamics.stimStart = interp1(t,1:length(t),0,'nearest'); % Find zero of data which reflects some task condition we made
rawreactionTime = Behaviour.reactionTime;
neuralDynamics.mreactionTime = interp1(t,1:length(t),mean(rawreactionTime),'nearest');
reactionTime = interp1(t,1:length(t),rawreactionTime,'nearest');
reactionTime(isnan(reactionTime)) = length(t);

PQ = [];
idx = discretize(reactionTime,20);

for n = 1:length(reactionTime)
    dat = meanTraj(X,n,6); %grab and calculate distance per reaction time
    x = dat(1,:); y = dat(2,:); z = dat(3,:);
    neuralDynamics.PQ(n) = sqrt((x(end)-x(1))^2+(y(end)-y(1))^2+(z(end)-z(1))^2)/rawreactionTime(n); %Calculate Distance
end
neuralDynamics.rawreactionTime = rawreactionTime;

%% Grab wave data
waveDynamics = struct();
if ~isempty(Waves)
    if ~isfield(Waves.wavesHit, 's')
        waveDynamics = getWaveDynamics2(Waves,Behaviour.parameters);
    else
        waveDynamics = getWaveDynamics(Waves);
    end
    
end

% %% Statistics of Neural Traj and Waves
% X = neuralTrajHitMiss;
% PQ = [];
% idx = discretize(reactionTime,20);
% 
% for n = 1:length(hittrials)
%     dat = meanTraj(X,n,6); %grab and calculate distance per reaction time
%     x = dat(1,1:reactionTime(n)); y = dat(2,1:reactionTime(n)); z = dat(3,1:reactionTime(n));
%     PQ(n) = sqrt((x(end)-x(1))^2+(y(end)-y(1))^2+(z(end)-z(1))^2); %Calculate Speed
%     dat1 = vertcat(zeros(1,6),rprimeh)'; % Grab norm vector alignment data
%     x = dat1(1,1:reactionTime(n)); y = dat1(2,1:reactionTime(n)); z = dat1(3,1:reactionTime(n));
%     TrajAngle(n) = mean(mean([x;y;z])); % Calculates how well aligned trajectories are for all trials
%     if ~isempty(Waves)
%         dPGD(n) = mean(diff(wavePGDhit(n,stimStart:reactionTime(n))));
%     end
% end
% %%
% X = neuralTrajHitMiss;
% PQ = [];
% idx = discretize(reactionTime,20);
% 
% for n = 1:length(hittrials)
%     dat = meanTraj(X,n,6); %grab and calculate distance per reaction time
%     x = dat(1,:); y = dat(2,:); z = dat(3,:);
%     PQ(n) = sqrt((x(end)-x(1))^2+(y(end)-y(1))^2+(z(end)-z(1))^2); %Calculate Distance
%     dat1 = vertcat(zeros(1,6),rprimeh)'; % Grab norm vector alignment data
%     x = dat1(1,1:reactionTime(n)); y = dat1(2,1:reactionTime(n)); z = dat1(3,1:reactionTime(n));
%     TrajAngle(n) = mean(mean([x;y;z])); % Calculates how well aligned trajectories are for all trials
%     if ~isempty(Waves)
%         dPGD(n) = mean(diff(wavePGDhit(n,stimStart:reactionTime(n))));
%     end
% end
% %%
% PQ = (PQ./rawreactionTime)';
% 
% figure,scatter(TrajAngle,rawreactionTime,'k','filled')
% xlim([0 .02])
% ylim([0 2])
% set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
% xlabel('Neural trajectory alignment'),ylabel('Reaction time (s)')
% 
% figure,scatter(PQ,rawreactionTime,'k','filled')
% xlim([0 15])
% ylim([0 2])
% set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
% xlabel('Neural trajectory speed'),ylabel('Reaction time (s)')
% 
% mdl = fitlm(PQ,rawreactionTime)
% figure,plot(mdl)
% xlim([0 15])
% ylim([0 2])
% set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
% xlabel('Neural trajectory speed'),ylabel('Reaction time (s)')
% legend(num2str(mdl.Rsquared.Ordinary),num2str(mdl.Coefficients.pValue(2)))
% title('')
% %%
% figure,scatter(dPGD,TrajAngle,'k','filled')
% set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
% xlabel('Neural trajectory speed'),ylabel('Phase gradient')
% xlim([0 10])
% 
% figure,scatter(PQ,dPGD,'k','filled')
% set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
% xlabel('Neural trajectory speed'),ylabel('Phase gradient')
% xlim([0 10])
% 
% mdl = fitlm(PQ,dPGD)
% figure,plot(mdl)
% set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
% xlabel('Neural trajectory speed'),ylabel('Phase gradient')
% legend(num2str(mdl.Rsquared.Ordinary),num2str(mdl.Coefficients.pValue(2)))
% xlim([0 10])
% title('')
% 
% figure,scatter(dPGD,rawreactionTime,'k','filled')
% set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
% xlabel('Phase gradient'),ylabel('Reaction time (s)')
% 
% mdl = fitlm(dPGD,rawreactionTime)
% figure,plot(mdl)
% set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
% xlabel('Phase gradient'),ylabel('Reaction time (s)')
% legend(num2str(mdl.Rsquared.Ordinary),num2str(mdl.Coefficients.pValue(2)))

end

%% Local functions

function [r,zscore_s,stability,matrix] = meanTraj(X,trials,components)
matrix = X(1:components,:,trials);
r = squeeze(mean(matrix,3));
% Subtract the mean from each element
centered_matrix = matrix - r;

% Calculate the Euclidean distance for each row
s = squeeze(sqrt(sum(centered_matrix.^2, 1)));
end

function speed_struct = speedTraj(X,trials,components,behaviour)
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

% Define behavioral state indices
pre_cue_idx = 1:75;
pre_movement_idx = 76:(76+ceil(mean(behaviour.reactionTime)*20));
reward_idx = (76+ceil(mean(behaviour.reactionTime)*20)):150;

% Calculate average speed for each state across all trials
% avg_speed_pre_cue = mean(speed_data(:,pre_cue_idx, :), [2 3]);
% avg_speed_pre_movement = mean(speed_data(:,pre_movement_idx, :), [2 3]);
% avg_speed_reward = mean(speed_data(:,reward_idx, :), [2 3]);

speed_struct.speed = speed_data;
speed_struct.preCue = speed_data(:,pre_cue_idx, :);
speed_struct.preMovement = speed_data(:,pre_movement_idx, :);
speed_struct.reward = speed_data(:,reward_idx, :);
end

function [neuralTrajSim,neuralTrajdiff,rprimehnorm,rprimemnorm] = neuralTrajDiff(r1,r2,varargin)
if strcmp(varargin,'initial'),initCalc = 1;else, initCalc = 0;end
if initCalc %r'(t)/||r'(t)|| * r(0)-ri(t)/||r(0)-ri(t)||
    rprimeh = diff(r1);
    rprimehnorm = rprimeh/norm(rprimeh);
    rprimem = r2(1,:)-r1; % control condition for reference
    rprimemnorm = rprimem/norm(rprimem);
    neuralTrajSim = dot(rprimehnorm',rprimemnorm(1:end-1,:)');
    neuralTrajdiff = rprimeh-rprimem;
else
    rprimeh = diff(r1);
    rprimehnorm = rprimeh/norm(rprimeh);
    rprimem = diff(r2);
    rprimemnorm = rprimem/norm(rprimem);
    neuralTrajSim = dot(rprimehnorm',rprimemnorm');
    neuralTrajdiff = rprimeh-rprimem;
end
end

function waveDynamics = getWaveDynamics(Waves)
rawWaveDensityhit = arrayfun(@(x) vertcat(x.wavePresent),Waves.wavesHit,'UniformOutput',false);waveDynamics.rawWaveDensityhit= vertcat(rawWaveDensityhit{:});
rawWavePGDhit = arrayfun(@(x) vertcat(x.PGD), Waves.wavesHit,'UniformOutput',false);waveDynamics.rawWavePGDhit = vertcat(rawWavePGDhit{:});
rawWaveSpeedhit = arrayfun(@(x) vertcat(x.s), Waves.wavesHit,'UniformOutput',false);waveDynamics.rawWaveSpeedhit = vertcat(rawWaveSpeedhit{:});

rawWaveDensitymiss = arrayfun(@(x) vertcat(x.wavePresent),Waves.wavesMiss,'UniformOutput',false);waveDynamics.rawWaveDensitymiss= vertcat(rawWaveDensitymiss{:});
rawWavePGDmiss = arrayfun(@(x) vertcat(x.PGD), Waves.wavesMiss,'UniformOutput',false);waveDynamics.rawWavePGDmiss = vertcat(rawWavePGDmiss{:});
rawWaveSpeedmiss = arrayfun(@(x) vertcat(x.s), Waves.wavesMiss,'UniformOutput',false);waveDynamics.rawWaveSpeedmiss = vertcat(rawWaveSpeedmiss{:});

rawWaveDensityFA = arrayfun(@(x) vertcat(x.wavePresent),Waves.wavesMIFA,'UniformOutput',false);waveDynamics.rawWaveDensityFA= vertcat(rawWaveDensityFA{:});
rawWavePGDFA = arrayfun(@(x) vertcat(x.PGD), Waves.wavesMIFA,'UniformOutput',false);waveDynamics.rawWavePGDFA = vertcat(rawWavePGDFA{:});
rawWaveSpeedFA = arrayfun(@(x) vertcat(x.s), Waves.wavesMIFA,'UniformOutput',false);waveDynamics.rawWaveSpeedFA = vertcat(rawWaveSpeedFA{:});

rawWaveDensityMIHit = arrayfun(@(x) vertcat(x.wavePresent),Waves.wavesMIHit,'UniformOutput',false);waveDynamics.rawWaveDensityMIHit= vertcat(rawWaveDensityMIHit{:});
rawWavePGDMIHit= arrayfun(@(x) vertcat(x.PGD), Waves.wavesMIHit,'UniformOutput',false);waveDynamics.rawWavePGDMIHit = vertcat(rawWavePGDMIHit{:});
rawWaveSpeedMIHit = arrayfun(@(x) vertcat(x.s), Waves.wavesMIHit,'UniformOutput',false);waveDynamics.rawWaveSpeedMIHit = vertcat(rawWaveSpeedMIHit{:});


% waveDensityhit = [];
% wavePGDhit = [];
% waveSpeedhit = [];
% 
% waveDensitymiss = [];
% wavePGDmiss = [];
% waveSpeedmiss = [];
% 
% waveDensityFA = [];
% wavePGDFA = [];
% waveSpeedFA = [];
% 
% waveDensityMIHit = [];
% wavePGDMIHit = [];
% waveSpeedMIHit = [];

% win = ceil(1:20:size(rawWaveDensityhit,2));
% for n = 1:length(win)-1
%     waveDensityhit = horzcat(waveDensityhit,sum(rawWaveDensityhit(:,win(n):win(n+1)),2)); %convert to wave/sec
%     wavePGDhit = horzcat(wavePGDhit,mean(rawWavePGDhit(:,win(n):win(n+1)),2));
%     waveSpeedhit = horzcat(waveSpeedhit,mean(rawWaveSpeedhit(:,win(n):win(n+1)),2));
%     
%     waveDensitymiss = horzcat(waveDensitymiss,sum(rawWaveDensitymiss(:,win(n):win(n+1)),2)); %convert to wave/sec
%     wavePGDmiss = horzcat(wavePGDmiss,mean(rawWavePGDmiss(:,win(n):win(n+1)),2));
%     waveSpeedmiss = horzcat(waveSpeedmiss,mean(rawWaveSpeedmiss(:,win(n):win(n+1)),2));
%     
%     waveDensityFA = horzcat(waveDensityFA,sum(rawWaveDensityFA(:,win(n):win(n+1)),2)); %convert to wave/sec
%     wavePGDFA = horzcat(wavePGDFA,mean(rawWavePGDFA(:,win(n):win(n+1)),2));
%     waveSpeedFA = horzcat(waveSpeedFA,mean(rawWaveSpeedFA(:,win(n):win(n+1)),2));
%     
%     waveDensityMIHit = horzcat(waveDensityMIHit,sum(rawWaveDensityMIHit(:,win(n):win(n+1)),2)); %convert to wave/sec
%     wavePGDMIHit = horzcat(wavePGDMIHit,mean(rawWavePGDMIHit(:,win(n):win(n+1)),2));
%     waveSpeedMIHit = horzcat(waveSpeedMIHit,mean(rawWaveSpeedMIHit(:,win(n):win(n+1)),2));
% end
end



function waveDynamics = getWaveDynamics2(Waves,parameters)

nPoints = 30; interval = (parameters.Fs*(parameters.windowAfterCue+parameters.windowBeforeCue))/nPoints;
for i=1:nPoints
    st = (i-1)*interval + 1;
    sp = (i)*interval + 1;
    waveDynamics.rawWaveSpeedHit(i).speed = horzcat(selectWaves(Waves.wavesHit,st,sp).speed);
    waveDynamics.rawWaveSpeedMiss(i).speed = horzcat(selectWaves(Waves.wavesMiss,st,sp).speed);
    waveDynamics.rawWaveSpeedMIHit(i).speed = horzcat(selectWaves(Waves.wavesMIHit,st,sp).speed);
    waveDynamics.rawWaveSpeedMIFA(i).speed = horzcat(selectWaves(Waves.wavesMIFA,st,sp).speed);
end

rawWaveDensityhit = arrayfun(@(x) vertcat(x.wavePresent),Waves.wavesHit,'UniformOutput',false);waveDynamics.rawWaveDensityhit= vertcat(rawWaveDensityhit{:});
rawWavePGDhit = arrayfun(@(x) vertcat(x.PGD), Waves.wavesHit,'UniformOutput',false);waveDynamics.rawWavePGDhit = vertcat(rawWavePGDhit{:});

rawWaveDensitymiss = arrayfun(@(x) vertcat(x.wavePresent),Waves.wavesMiss,'UniformOutput',false);waveDynamics.rawWaveDensitymiss= vertcat(rawWaveDensitymiss{:});
rawWavePGDmiss = arrayfun(@(x) vertcat(x.PGD), Waves.wavesMiss,'UniformOutput',false);waveDynamics.rawWavePGDmiss = vertcat(rawWavePGDmiss{:});

rawWaveDensityFA = arrayfun(@(x) vertcat(x.wavePresent),Waves.wavesMIFA,'UniformOutput',false);waveDynamics.rawWaveDensityFA= vertcat(rawWaveDensityFA{:});
rawWavePGDFA = arrayfun(@(x) vertcat(x.PGD), Waves.wavesMIFA,'UniformOutput',false);waveDynamics.rawWavePGDFA = vertcat(rawWavePGDFA{:});

rawWaveDensityMIHit = arrayfun(@(x) vertcat(x.wavePresent),Waves.wavesMIHit,'UniformOutput',false);waveDynamics.rawWaveDensityMIHit= vertcat(rawWaveDensityMIHit{:});
rawWavePGDMIHit= arrayfun(@(x) vertcat(x.PGD), Waves.wavesMIHit,'UniformOutput',false);waveDynamics.rawWavePGDMIHit = vertcat(rawWavePGDMIHit{:});


% waveDensityhit = [];
% wavePGDhit = [];
% waveSpeedhit = [];
% 
% waveDensitymiss = [];
% wavePGDmiss = [];
% waveSpeedmiss = [];
% 
% waveDensityFA = [];
% wavePGDFA = [];
% waveSpeedFA = [];
% 
% waveDensityMIHit = [];
% wavePGDMIHit = [];
% waveSpeedMIHit = [];

% win = ceil(1:20:size(rawWaveDensityhit,2));
% for n = 1:length(win)-1
%     waveDensityhit = horzcat(waveDensityhit,sum(rawWaveDensityhit(:,win(n):win(n+1)),2)); %convert to wave/sec
%     wavePGDhit = horzcat(wavePGDhit,mean(rawWavePGDhit(:,win(n):win(n+1)),2));
%     waveSpeedhit = horzcat(waveSpeedhit,mean(rawWaveSpeedhit(:,win(n):win(n+1)),2));
%     
%     waveDensitymiss = horzcat(waveDensitymiss,sum(rawWaveDensitymiss(:,win(n):win(n+1)),2)); %convert to wave/sec
%     wavePGDmiss = horzcat(wavePGDmiss,mean(rawWavePGDmiss(:,win(n):win(n+1)),2));
%     waveSpeedmiss = horzcat(waveSpeedmiss,mean(rawWaveSpeedmiss(:,win(n):win(n+1)),2));
%     
%     waveDensityFA = horzcat(waveDensityFA,sum(rawWaveDensityFA(:,win(n):win(n+1)),2)); %convert to wave/sec
%     wavePGDFA = horzcat(wavePGDFA,mean(rawWavePGDFA(:,win(n):win(n+1)),2));
%     waveSpeedFA = horzcat(waveSpeedFA,mean(rawWaveSpeedFA(:,win(n):win(n+1)),2));
%     
%     waveDensityMIHit = horzcat(waveDensityMIHit,sum(rawWaveDensityMIHit(:,win(n):win(n+1)),2)); %convert to wave/sec
%     wavePGDMIHit = horzcat(wavePGDMIHit,mean(rawWavePGDMIHit(:,win(n):win(n+1)),2));
%     waveSpeedMIHit = horzcat(waveSpeedMIHit,mean(rawWaveSpeedMIHit(:,win(n):win(n+1)),2));
% end
end