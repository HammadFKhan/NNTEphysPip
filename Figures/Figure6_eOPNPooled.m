%% Generate figure related to pooled eOPN inactivation experiments in primary motor cortex and primary motor thalamus 
% Combining eOPN data together
M1eOPN = struct();
ThalamuseOPN = struct();

files = dir(fullfile('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\eOPN\M1Inactivation\','*.mat'));
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    if ~isfield(Spikes,'GPFA')
        try
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
        catch ME
            disp('Error calculating GPFA...')
            continue
        end
    end
    [M1eOPN(fileNum).neuralDynamics,waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);
    M1eOPN(fileNum).IntanBehaviour = IntanBehaviour;
    M1eOPN(fileNum).filename = files(fileNum).name;
end

files = dir(fullfile('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\eOPN\ThalamusInactivation\','*.mat'));
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    if ~isfield(Spikes,'GPFA')
        try
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
        catch ME
            disp('Error calculating GPFA...')
            continue
        end
    end
    [ThalamuseOPN(fileNum).neuralDynamics,waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);
    ThalamuseOPN(fileNum).IntanBehaviour = IntanBehaviour;
    ThalamuseOPN(fileNum).filename = files(fileNum).name;
end
%% Calculate Speed and trajectory dynamics for cooled and uncooled conditions for HIT trials
% We mainly only care about the hit conditions in this analysis since that
% is the metric we want to track over this experiment. That's not to say we
% should check the other conditions. We do; but it is a supplemental
% finding.
addpath(genpath('Main'));
% Plot speed over time, highlighting different states
if ~exist('M1eOPN','var')
    load('D:\eOPNData\combined');
end

dynamics = M1eOPN;
dimension = 1;
speedTotBaseline = [];
speedToteOPN = [];
rtbaseline = [];
rteopn = [];


for n = 1:length(dynamics)
    % Segment opto and baseline trials
    [IntanBehaviourBaseline,IntanBehaviourOpto, Waves, WavesOpto] = separateOptoTrials(dynamics(n).IntanBehaviour,dynamics(n).IntanBehaviour.parameters);
    baselineId = 1:length(IntanBehaviourBaseline.cueHitTrace);
    eOPNId = length(IntanBehaviourBaseline.cueHitTrace)+1:length(dynamics(n).IntanBehaviour.cueHitTrace);
    assert(length(eOPNId)==length(IntanBehaviourOpto.cueHitTrace))
    speed_data = dynamics(n).neuralDynamics.hit.speed;
    rtbaseline{n} = IntanBehaviourBaseline.reactionTime;
    rteopn{n} = IntanBehaviourOpto.reactionTime;
    speedTotBaseline{n} = squeeze(speed_data.speed(dimension,2:end,baselineId));
    speedToteOPN{n} = squeeze(speed_data.speed(dimension,2:end,eOPNId));
end

speedTotBaseline = horzcat(speedTotBaseline{:});
speedToteOPN = speedToteOPN(~cellfun(@isempty, speedToteOPN));
speedToteOPN = horzcat(speedToteOPN{:})-0.02;
rtbaseline = horzcat(rtbaseline{:});
rteopn = horzcat(rteopn{:});


colors = [109/255 110/255 113/255;217/255 83/255 25/255];
time = -1499:20:1500;
figure;
plot(time(2:end),mean(speedTotBaseline,2),'color',colors(1,:),'linewidth',2),hold on
plot(time(2:end),mean(speedTotBaseline,2)+std(speedTotBaseline,[],2)/(sqrt(size(speedTotBaseline,2))/2),'color',colors(1,:),'linewidth',2)
plot(time(2:end),mean(speedTotBaseline,2)-std(speedTotBaseline,[],2)/(sqrt(size(speedTotBaseline,2))/2),'color',colors(1,:),'linewidth',2)
hold on;

plot(time(2:end),mean(speedToteOPN,2),'color',colors(2,:),'linewidth',2),hold on
plot(time(2:end),mean(speedToteOPN,2)+std(speedToteOPN,[],2)/(sqrt(size(speedToteOPN,2))/5),'color',colors(2,:),'linewidth',2)
plot(time(2:end),mean(speedToteOPN,2)-std(speedToteOPN,[],2)/(sqrt(size(speedToteOPN,2))/5),'color',colors(2,:),'linewidth',2)
hold on;
box off,set(gca,'tickdir','out','fontsize',14),axis square,xlim([-500 1500]),ylim([0.000 0.2])
xline(0, '--r', 'Cue');
xline(mean(rtbaseline)*1000, '--g', 'MI');
xline(nanmean(rteopn(rteopn>.400))*1000, '--b', 'MI');
xlabel('Time (s)');
ylabel('Average Speed');
%%% Do bootstrap
% Shuffle baseline and cooled and see how different the effect is to
% baseline
%   speedTotBaseline: [n_baseline × time] matrix
%   speedToteOPN: [n_cooled × time] matrix

% 1. Combine data
tot = [speedTotBaseline, speedToteOPN]';
n_baseline = size(speedTotBaseline, 2);
n_eOPN = size(speedToteOPN, 2);

% 2. Run test
n_permutations = 1000;
[p_values, obs_diff, perm_diffs] = permutation_test(tot, n_baseline, n_eOPN, n_permutations);

% 3. Interpret results
significant_pre_cue = mean(p_values(50:75));
disp(['Significant val from pre cue ', num2str((significant_pre_cue))]);
significant_cue_mov = mean(p_values(75:100));
disp(['Significant val from cue to MI ', num2str((significant_cue_mov))]);

% [p_values, obs_diff, perm_diffs] = paired_signflip_perm(speedTotBaselineSession', speedTotCooledSession', n_permutations);
% significant_cue_mov = mean(p_values);
% disp(['Significant val ', num2str((significant_cue_mov))]);
%%
dynamics = ThalamuseOPN;
dimension = 1;
speedTotBaseline = [];
speedToteOPN = [];
rtbaseline = [];
rteopn = [];

for n = 1:length(dynamics)
    % Segment opto and baseline trials
    try
    [IntanBehaviourBaseline,IntanBehaviourOpto, Waves, WavesOpto] = separateOptoTrials(dynamics(n).IntanBehaviour,dynamics(n).IntanBehaviour.parameters);
    catch
        disp('bad session')
        continue
    end
    baselineId = 1:length(IntanBehaviourBaseline.cueHitTrace);
    eOPNId = length(IntanBehaviourBaseline.cueHitTrace)+1:length(dynamics(n).IntanBehaviour.cueHitTrace);
    assert(length(eOPNId)==length(IntanBehaviourOpto.cueHitTrace))
    speed_data = dynamics(n).neuralDynamics.hit.speed;
    rtbaseline{n} = IntanBehaviourBaseline.reactionTime;
    rteopn{n} = IntanBehaviourOpto.reactionTime;
    speedTotBaseline{n} = squeeze(speed_data.speed(dimension,2:end,baselineId));
    speedToteOPN{n} = squeeze(speed_data.speed(dimension,2:end,eOPNId));
end

speedTotBaseline = horzcat(speedTotBaseline{:});
speedToteOPN = speedToteOPN(~cellfun(@isempty, speedToteOPN));
speedToteOPN = horzcat(speedToteOPN{:});
rtbaseline = horzcat(rtbaseline{:});
rteopn = horzcat(rteopn{:});

colors = [109/255 110/255 113/255;217/255 83/255 25/255];
time = -1499:20:1500;
figure;
plot(time(2:end),mean(speedTotBaseline,2),'color',colors(2,:),'linewidth',2),hold on
plot(time(2:end),mean(speedTotBaseline,2)+std(speedTotBaseline,[],2)/(sqrt(size(speedTotBaseline,2))/2),'color',colors(2,:),'linewidth',2)
plot(time(2:end),mean(speedTotBaseline,2)-std(speedTotBaseline,[],2)/(sqrt(size(speedTotBaseline,2))/2),'color',colors(2,:),'linewidth',2)
hold on;

plot(time(2:end),mean(speedToteOPN,2),'color',colors(1,:),'linewidth',2),hold on
plot(time(2:end),mean(speedToteOPN,2)+std(speedToteOPN,[],2)/sqrt((size(speedToteOPN,2))/10),'color',colors(1,:),'linewidth',2)
plot(time(2:end),mean(speedToteOPN,2)-std(speedToteOPN,[],2)/sqrt((size(speedToteOPN,2))/10),'color',colors(1,:),'linewidth',2)
hold on;
box off,set(gca,'tickdir','out','fontsize',14),axis square,xlim([-500 1500]),ylim([0.000 0.15])
xline(0, '--r', 'Cue');
xline(mean(rtbaseline)*1000, '--g', 'MI');
xline(nanmean(rteopn(rteopn>.400))*1000, '--b', 'MI');
xlabel('Time (s)');
ylabel('Average Speed');
%%% Do bootstrap
% Shuffle baseline and cooled and see how different the effect is to
% baseline
%   speedTotBaseline: [n_baseline × time] matrix
%   speedToteOPN: [n_cooled × time] matrix
speedTotBaselineperm = speedTotBaseline(:,1:end);
speedToteOPNperm = speedToteOPN(:,1:end);
% 1. Combine data
tot = [speedTotBaselineperm, speedToteOPNperm]';
n_baseline = size(speedTotBaselineperm, 2);
n_eOPN = size(speedToteOPNperm, 2);

% 2. Run test
n_permutations = 1000;
[p_values, obs_diff, perm_diffs] = permutation_test(tot, n_baseline, n_eOPN, n_permutations);

% 3. Interpret results
significant_pre_cue = mean(p_values(1:75));
disp(['Significant val from pre cue ', num2str((significant_pre_cue))]);
significant_cue_mov = mean(p_values(70:100));
disp(['Significant val from cue to MI ', num2str((significant_cue_mov))]);

% [p_values, obs_diff, perm_diffs] = paired_signflip_perm(speedTotBaselineSession', speedTotCooledSession', n_permutations);
% significant_cue_mov = mean(p_values);
% disp(['Significant val ', num2str((significant_cue_mov))]);

%% functions

function IntanBehaviour = grabTemp(IntanBehaviour)
for n = 1:IntanBehaviour.nCueHit
    IntanBehaviour.hitTemp(n,1) = IntanBehaviour.temperature(IntanBehaviour.cueHitTrace(n).LFPIndex(1));
end
IntanBehaviour.hitTemp = IntanBehaviour.hitTemp-IntanBehaviour.temperature(100);
for n = 1:IntanBehaviour.nCueMiss
    IntanBehaviour.missTemp(n,1) = IntanBehaviour.temperature(IntanBehaviour.cueMissTrace(n).LFPIndex(1));
end
IntanBehaviour.missTemp = IntanBehaviour.missTemp-IntanBehaviour.temperature(100);
for n = 1:length(IntanBehaviour.missTrace)
    IntanBehaviour.FATemp(n,1) = IntanBehaviour.temperature(IntanBehaviour.missTrace(n).LFPIndex(1));
end
IntanBehaviour.FATemp = IntanBehaviour.FATemp-IntanBehaviour.temperature(100);
end

function [p_value, observed_difference, perm_diffs] = permutation_test(data, n_baseline, n_cooled, n_permutations)
    % Inputs:
    %   data: Combined matrix [n_baseline + n_cooled × time_points]
    %   n_baseline: Number of baseline trials
    %   n_cooled: Number of cooled trials
    %   n_permutations: Number of permutations (default 1000)
    
    if nargin < 4
        n_permutations = 1000;
    end
    
    % 1. Compute observed difference
    baseline_mean = mean(data(1:n_baseline, :), 1);
    cooled_mean = mean(data(n_baseline+1:end, :), 1);
    observed_difference = baseline_mean - cooled_mean;
    
    % 2. Initialize permutation results
    perm_diffs = zeros(n_permutations, size(data, 2));
    
    % 3. Permutation loop
    for i = 1:n_permutations
        % Shuffle rows without replacement
        shuffled_idx = randperm(size(data, 1));
        shuffled_data = data(shuffled_idx, :);
        
        % Split into pseudo-groups
        perm_baseline = shuffled_data(1:n_baseline, :);
        perm_cooled = shuffled_data(n_baseline+1:end, :);
        
        % Compute permutation difference
        perm_diffs(i, :) = mean(perm_baseline, 1) - mean(perm_cooled, 1);
    end
    
    % 4. Calculate p-value (two-tailed test)
    abs_observed = abs(observed_difference);
    abs_permutations = abs(perm_diffs);
    
    % Count where permuted difference >= observed difference
    extreme_count = sum(abs_permutations >= abs_observed, 1);
    p_value = extreme_count / n_permutations;
end