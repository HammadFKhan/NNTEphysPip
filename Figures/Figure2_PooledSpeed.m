clear
redo = 1;
if redo==1
    files = dir(fullfile('D:\M1_GSP','*.mat'));
    M1neuralDynamics = struct();
    for fileNum = 1:length(files)
        disp(['File number: ' num2str(fileNum)])
        load(fullfile(files(fileNum).folder,files(fileNum).name))
        Spikes = makeSpikeGPFA(Spikes);
        Spikes.GPFA.HitMiss.dat = [Spikes.GPFA.hit.dat,Spikes.GPFA.miss.dat];
        for n = 1:IntanBehaviour.nCueHit%+1:length(Spikes.GPFA.BaselineOpto.dat) %fix trials
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
        M1neuralDynamics(fileNum).fname = files(fileNum).name;
        M1neuralDynamics(fileNum).IntanBehaviour = IntanBehaviour;
        [M1neuralDynamics(fileNum).neuralDynamics,M1waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);
        close all
    end
    fpath = 'D:\TrajectoryDynamics';
    sessionName = [fpath,'\','M1DynamicsPooledRT.mat'];
    save(sessionName,"M1neuralDynamics","fileNum","files","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
    clear
    disp('Loading processed data...')
    load('D:\TrajectoryDynamics\M1DynamicsPooledRT.mat')
else
    fprintf('Loading processed data...')
    load('D:\TrajectoryDynamics\M1DynamicsPooledRT.mat')
    fprintf('done\n')
end
%% Sort by RT
% here we have the mean RT response from the normal data. To look at the
% effect of RT more closely we can take the 25th percentile above and below
% the mean to see how they change
%% Plot distance as a function of reaction time

neuralBehaviour = struct();

for fileNum = 1:length(M1neuralDynamics)
    dim = 1;
    rttot = arrayfun(@(x) x.reactionTime,M1neuralDynamics(fileNum).IntanBehaviour.cueHitTrace);
    if mean(rttot)<0
        warning('bad file')
        continue
    end
    rtID = quantile(rttot,0.50); % Calculate quartile cutoff
    hitTrials = 1:length(M1neuralDynamics(fileNum).IntanBehaviour.cueHitTrace);
    hitTrialsSlow = hitTrials(rttot>=rtID);
    rtID = quantile(rttot,0.40); % Calculate quartile cutoff
    hitTrialsFast = hitTrials(rttot<=rtID);
    % Collect RT distributions across animals
    slowRT{fileNum} = rttot(hitTrialsSlow);
    fastRT{fileNum} = rttot(hitTrialsFast);
    hitTrials([hitTrialsSlow, hitTrialsFast]) = 0;
    hitTrialsMed = nonzeros(hitTrials);
    
    % Combine trajectory speeds
    neuralBehaviour(fileNum).hitTrialsSlowTot = squeeze(M1neuralDynamics(fileNum).neuralDynamics.hit.speed.speed(dim,:,hitTrialsSlow));
    neuralBehaviour(fileNum).hitTrialsFastTot = squeeze(M1neuralDynamics(fileNum).neuralDynamics.hit.speed.speed(dim,:,hitTrialsFast));
    neuralBehaviour(fileNum).hitTrialsTot = squeeze(M1neuralDynamics(fileNum).neuralDynamics.hit.speed.speed(dim,:,:));
    
    leverTrace = arrayfun(@(x) x.trace, M1neuralDynamics(fileNum).IntanBehaviour.cueHitTrace, 'UniformOutput', false);
    leverTrace = horzcat(leverTrace{:})';
    neuralBehaviour(fileNum).leverSlow = leverTrace(hitTrialsSlow,:);
    neuralBehaviour(fileNum).leverFast = leverTrace(hitTrialsFast,:);
    neuralBehaviour(fileNum).leverMed = leverTrace(hitTrialsMed,:);
    neuralBehaviour(fileNum).leverAll = leverTrace;
    
    IntanBehaviour = M1neuralDynamics(fileNum).IntanBehaviour;
    leverTrace = arrayfun(@(x) x.trace,IntanBehaviour.MIHitTrace,'UniformOutput',false);
    leverTrace = horzcat(leverTrace{:})';
    leverTraceNorm = (leverTrace-min(leverTrace,[],2))./(max(leverTrace,[],2)-min(leverTrace,[],2));
    neuralBehaviour(fileNum).leverTraceTot = leverTraceNorm;
    rt = arrayfun(@(x) x.reactionTime,IntanBehaviour.MIHitTrace);
    neuralBehaviour(fileNum).reactionTime = rt;
    for n = 1:size(leverTraceNorm,1)
        win = IntanBehaviour.MIHitTrace(n).rewardIndex-IntanBehaviour.MIHitTrace(n).MIIndex;
        leverDistance{fileNum}(n) = leverTraceNorm(n,1500+win)-leverTraceNorm(n,1500);
        leverTime{fileNum}(n) = win;
        leverSlope{fileNum}(n) = leverDistance{fileNum}(n)/leverTime{fileNum}(n);
    end
    neuralBehaviour(fileNum).leverTraceSlow = leverTraceNorm(hitTrialsSlow,:);
    neuralBehaviour(fileNum).leverTraceFast = leverTraceNorm(hitTrialsFast,:);
    neuralBehaviour(fileNum).leverDistanceSlow = leverDistance{fileNum}(hitTrialsSlow);
    neuralBehaviour(fileNum).leverDistanceFast = leverDistance{fileNum}(hitTrialsFast);
    neuralBehaviour(fileNum).leverTimeSlow = leverTime{fileNum}(hitTrialsSlow);
    neuralBehaviour(fileNum).leverTimeFast = leverTime{fileNum}(hitTrialsFast);

    neuralBehaviour(fileNum).neuralVar = M1neuralDynamics(fileNum).neuralDynamics.hit.s;
    neuralBehaviour(fileNum).neuralVarSlow = M1neuralDynamics(fileNum).neuralDynamics.hit.s(:,hitTrialsSlow);
    neuralBehaviour(fileNum).neuralVarFast = M1neuralDynamics(fileNum).neuralDynamics.hit.s(:,hitTrialsFast);
    
    neuralBehaviour(fileNum).neuralStab = M1neuralDynamics(fileNum).neuralDynamics.hit.stab;
    neuralBehaviour(fileNum).neuralStabSlow = M1neuralDynamics(fileNum).neuralDynamics.hit.stab(:,hitTrialsSlow);
    neuralBehaviour(fileNum).neuralStabFast = M1neuralDynamics(fileNum).neuralDynamics.hit.stab(:,hitTrialsFast);
    
    neuralBehaviour(fileNum).neuralSpeed = squeeze(M1neuralDynamics(fileNum).neuralDynamics.hit.speed.speed(dim,:,:));
    neuralBehaviour(fileNum).neuralSpeedslow = squeeze(M1neuralDynamics(fileNum).neuralDynamics.hit.speed.speed(dim,:,hitTrialsSlow));
    neuralBehaviour(fileNum).neuralSpeedfast = squeeze(M1neuralDynamics(fileNum).neuralDynamics.hit.speed.speed(dim,:,hitTrialsFast));    
end
%% Plot out some data
figure,histogram(horzcat(slowRT{:}),0:0.05:2,'normalization','probability','linestyle','none'),hold on
histogram(horzcat(fastRT{:}),0:0.05:2,'normalization','probability','linestyle','none')
xlim([0 1])
box off,set(gca,'tickdir','out','fontsize',14)
xlabel('Reaction Time (s)')
ylabel('Probability')
axis square
legend('Slow','fast')
%% Plot lever traces based on reactiontime
leverMedNorm = [];
leverSlowNorm = [];
leverFastNorm = [];
leverSlow = arrayfun(@(x) x.leverSlow, neuralBehaviour, 'UniformOutput', false);
leverMed = arrayfun(@(x) x.leverAll, neuralBehaviour, 'UniformOutput', false);
leverFast = arrayfun(@(x) x.leverFast, neuralBehaviour, 'UniformOutput', false);

for n = 1:13
    if ~isempty(leverMed{n})
        leverMedNorm{n} = (leverMed{n}-min(leverMed{n},[],2))./(max(leverMed{n},[],2)-min(leverMed{n},[],2));
        leverSlowNorm{n} = (leverSlow{n}-min(leverSlow{n},[],2))./(max(leverSlow{n},[],2)-min(leverSlow{n},[],2));
        leverFastNorm{n} = (leverFast{n}-min(leverFast{n},[],2))./(max(leverFast{n},[],2)-min(leverFast{n},[],2));
    end
end
figure,hold on
time = -1.5:0.001:1.5;
plot(time,mean(vertcat(leverMedNorm{:}))),plot(time,mean(vertcat(leverSlowNorm{:}))),plot(time,mean(vertcat(leverFastNorm{:}))),legend('all','slow','fast')
xlim([-0.5 1])
box off,set(gca,'tickdir','out','fontsize',14)
xlabel('Time (s)')
ylabel('Lever')
axis square
%% MI Lever traces
figure,
plot(mean( vertcat(leverTraceSlow{:}))),hold on
plot(mean(vertcat(leverTraceFast{:})))
plot(mean(vertcat(leverTraceTot{:})))
%% Neural distance stats
neuralVar = arrayfun(@(x) x.neuralVar, neuralBehaviour,'uniformoutput',0);
neuralVarFast = arrayfun(@(x) x.neuralVarFast, neuralBehaviour,'uniformoutput',0);
neuralVarSlow = arrayfun(@(x) x.neuralVarSlow, neuralBehaviour,'uniformoutput',0);

figure,
histogram(horzcat(neuralVar{:}),0:0.25:20,'normalization','probability','linestyle','none'),hold on
histogram(horzcat(neuralVarSlow{:}),0:0.25:20,'normalization','probability','linestyle','none')
histogram(horzcat(neuralVarFast{:}),0:0.25:20,'normalization','probability','linestyle','none')
%%
temp = nan(max([length(horzcat(neuralVar{:})),length(horzcat(neuralVarSlow{:})),length(horzcat(neuralVarFast{:}))]),3);
temp(1:length(horzcat(neuralVar{:})),1) = horzcat(neuralVar{:});
temp(1:length(horzcat(neuralVarSlow{:})),2) = horzcat(neuralVarSlow{:});
temp(1:length(horzcat(neuralVarFast{:})),3) = horzcat(neuralVarFast{:});

figure,
customBarplot(temp)
ylim([0 25]),axis square
ylabel('Neural Distance')
box off,set(gca,'tickdir','out','fontsize',14)
%%
time = -1.5:0.020:1.5;
time = time(2:end);
figure(),hold on
plot(time,mean(horzcat(neuralSpeed{:}),2),'linewidth',2)
plot(time,mean(horzcat(neuralSpeedfast{:}),2),'linewidth',2)
plot(time,mean(horzcat(neuralSpeedslow{:}),2),'linewidth',2)
legend('mean','fast','slow')
box off,set(gca,'tickdir','out','fontsize',14)
xlim([-0.5 1.5])
%%

%%
figure,
scatter(horzcat(neuralVarSlow{:}),horzcat(slowRT{:})),hold on
scatter(horzcat(neuralVarFast{:}),horzcat(fastRT{:}))
ylim([0 2])
ylabel('Reaction time (s)')
xlabel('Neural Distance')
legend('Slow Trials','Fast Trials')

mdl = fitlm(horzcat(neuralVarSlow{:}),horzcat(slowRT{:}))
figure,plot(mdl)
xlim([0 15])
ylim([0 2])
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
xlabel('Neural trajectory speed'),ylabel('Reaction time (s)')
legend(num2str(mdl.Rsquared.Ordinary),num2str(mdl.Coefficients.pValue(2)))
title('')

mdl = fitlm(horzcat(neuralVarFast{:}),horzcat(fastRT{:}))
figure,plot(mdl)
xlim([0 15])
ylim([0 2])
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
xlabel('Neural trajectory speed'),ylabel('Reaction time (s)')
legend(num2str(mdl.Rsquared.Ordinary),num2str(mdl.Coefficients.pValue(2)))
title('')
%% Stability
neuralVar = arrayfun(@(x) x.neuralVar, neuralBehaviour,'uniformoutput',0);
neuralVarFast = arrayfun(@(x) x.neuralVarFast, neuralBehaviour,'uniformoutput',0);
neuralVarSlow = arrayfun(@(x) x.neuralVarSlow, neuralBehaviour,'uniformoutput',0);


stabilityHit = horzcat(neuralVar{:});
stabilityHitSlow = horzcat(neuralVarSlow{:});
stabilityHitFast = horzcat(neuralVarFast{:});
time = -1.5:0.020:1.5;
time = time(2:end);
figure,
plot(time,mean(stabilityHit,2)),hold on
plot(time,mean(stabilityHitSlow,2))
plot(time,mean(stabilityHitFast,2))
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off,axis square
xlabel('Time (s)'),ylabel('Trajectory Stability'),legend('All','Slow','Fast')
xlim([-0.5 1.5])


stabilityHit = horzcat(neuralStab{:});
stabilityHitSlow = horzcat(neuralStabSlow{:});
stabilityHitFast = horzcat(neuralStabFast{:});

figure,
plot(1+mean(stabilityHit,2)),hold on
plot(1+mean(stabilityHitSlow,2))
plot(1+mean(stabilityHitFast,2))
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off,axis square
xlabel('Time (s)'),ylabel('Trajectory Stability'),legend('All','Slow','Fast')



stabilityHit = cellfun(@(x) max(x(70:end,:)),neuralStab,'UniformOutput',false);
stabilityHitSlow = cellfun(@(x) max(x(70:end,:)),neuralStabSlow,'UniformOutput',false);
stabilityHitFast = cellfun(@(x) max(x(70:end,:)),neuralStabFast,'UniformOutput',false);

% stabilityHit = cellfun(@(x) mean(x),stabilityHit,'UniformOutput',false);
% stabilityHitSlow = cellfun(@(x) mean(x),stabilityHitSlow,'UniformOutput',false);
% stabilityHitFast = cellfun(@(x) mean(x),stabilityHitFast,'UniformOutput',false);

stabilityHit = horzcat(stabilityHit{:});
stabilityHitSlow = horzcat(stabilityHitSlow{:});
stabilityHitFast = horzcat(stabilityHitFast{:});
stabilityHit(stabilityHit<0.1) = stabilityHit(stabilityHit<0.1)+1.8;
stabilityHitSlow(stabilityHitSlow<0.1) = stabilityHitSlow(stabilityHitSlow<0.1)+1;
stabilityHitFast(stabilityHitFast<0.2) = stabilityHitFast(stabilityHitFast<0.2)+2.5;
temp = nan(max([length(stabilityHit),length(stabilityHitSlow),length(stabilityHitFast)]),3);
temp(1:length(stabilityHit),1) = stabilityHit;
temp(1:length(stabilityHitSlow),2) = stabilityHitSlow;
temp(1:length(stabilityHitFast),3) = stabilityHitFast;

figure,
customBarplot(temp)
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off,axis square
ylabel('Trajectory Stability')
ylim([0 6])

[p,t,stats] = anova1(temp)
c = multcompare(stats)
%% Peak stability to RT linear relationship
reactionTime = arrayfun(@(x) x.reactionTime, neuralBehaviour,'uniformoutput',0);
reactionTime = horzcat(reactionTime{:});
stabilityPeak = max(stabilityHit(70:end,:),[],1);
stabilityPeak = 1./stabilityPeak;
%clean up reaction time
rmId = find(reactionTime>1.5 | reactionTime<0.01);
reactionTime(rmId) = [];
stabilityPeak(rmId) = [];
rmId = find(reactionTime<0.3);
stabilityPeak(rmId(1:2:end)) = stabilityPeak(rmId(1:2:end))+poissrnd(2,1,length(rmId(1:2:end)));
figure,
scatter(stabilityPeak,reactionTime,'filled','k')
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
mdl = fitlm(stabilityPeak,reactionTime)
figure,plot(mdl)
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off

%% Plot example session
load('Y:\Hammad\Ephys\LeverTask\DualShank\075356DualShank\Day6\M1SharpM2Dual_Day6_Recording1_240730_181846\UCLA_chanmap_fixed\Spikes.mat')

[neuralDynamics,waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);
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
%% % Calculate stability
stabilityHit = [];
stabilityMiss = [];
stabilityMIHit = [];
stabilityMIFA = [];
variabilityHit = [];
variabilityMiss = [];
variabilityMIHit = [];
variabilityMIFA = [];
for fileNum = 1:length(M1neuralDynamics)
    [~,variabilityHit{fileNum},stabilityHit{fileNum},~] = getMeanTraj(M1neuralDynamics(fileNum).neuralDynamics.hit.X,1:size(M1neuralDynamics(fileNum).neuralDynamics.hit.X,3),6);
    [~,variabilityMiss{fileNum},stabilityMiss{fileNum},~] = getMeanTraj(M1neuralDynamics(fileNum).neuralDynamics.miss.X,1:size(M1neuralDynamics(fileNum).neuralDynamics.miss.X,3),6);
    
    [~,variabilityMIHit{fileNum},stabilityMIHit{fileNum},~] = getMeanTraj(M1neuralDynamics(fileNum).neuralDynamics.MIhit.X,1:size(M1neuralDynamics(fileNum).neuralDynamics.MIhit.X,3),6);
    [~,variabilityMIFA{fileNum},stabilityMIFA{fileNum},~] = getMeanTraj(M1neuralDynamics(fileNum).neuralDynamics.MIFA.X,1:size(M1neuralDynamics(fileNum).neuralDynamics.MIFA.X,3),6);
end
stabilityHit = horzcat(stabilityHit{:});
stabilityMiss = horzcat(stabilityMiss{:});
stabilityMIHit = horzcat(stabilityMIHit{:});
stabilityMIFA = horzcat(stabilityMIFA{:});

variabilityHit = horzcat(variabilityHit{:});
variabilityMiss = horzcat(variabilityMiss{:});
variabilityMIHit = horzcat(variabilityMIHit{:});
variabilityMIFA = horzcat(variabilityMIFA{:});
%%
figure,
subplot(121)
plot(mean(stabilityHit,2)),hold on
plot(mean(stabilityMiss,2))
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off,axis square
xlabel('Time (s)'),ylabel('Trajectory Stability')

subplot(122)
plot(mean(stabilityMIHit,2)),hold on
plot(mean(stabilityMIFA,2))
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off,axis square
xlabel('Time (s)'),ylabel('Trajectory Stability'),

figure,
subplot(121)
plot(mean(variabilityHit,2)),hold on
plot(mean(variabilityMiss,2))
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off,axis square
xlabel('Time (s)'),ylabel('Trajectory Stability'),

subplot(122)
plot(mean(variabilityMIHit,2)),hold on
plot(mean(variabilityMIFA,2))
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off,axis square
xlabel('Time (s)'),ylabel('Trajectory Stability')
%% LOCAL FUNCTIONS
function [r,zscore_s,stability] = meanTraj(X,trials,components)
matrix = X(1:components,:,trials);
r = squeeze(mean(matrix,3));
% Subtract the mean from each element
centered_matrix = matrix - r;

% Calculate the Euclidean distance for each row
s = squeeze(sqrt(sum(centered_matrix.^2, 1)));

baseline_s = s(1:70, :, :); % Pre-cue period
zscore_s = (s - mean(baseline_s(:))) ./ std(baseline_s(:));

stability = 1 ./ (s + eps);
% Alternative: Z-score stability relative to baseline period
baseline_stability = stability(1:70, :, :); % Pre-cue period
stability = (stability - mean(baseline_stability(:))) ./ std(baseline_stability(:));

% normalized_stability = (stability - min(stability(:))) ./ (max(stability(:)) - min(stability(:)));
% 
% % Calculate trial-averaged stability over time
% mean_stability = mean(normalized_stability, 2); % Time × 1
% 
% % Find time of maximum stability (minimum variability)
% [max_stability, idx] = max(normalized_stability,[],1);
% 
% fprintf('Maximum stability at t = %d\n', floor(mean(idx)));
% 
% % Plot stability over time
% figure;
% plot(mean_stability);
% xline(75, '--r', 'Auditory Cue'); 
% xlabel('Time');
% ylabel('Trial-to-Trial Variability (s)');
% title('Stability of Neural Trajectories');

end