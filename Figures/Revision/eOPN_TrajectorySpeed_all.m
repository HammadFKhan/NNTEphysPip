
% M1 Dynamics for Th inactivation 
load('\\10.165.57.13\Sutter_backup\Hammad\Ephys\LeverTask\Data_for_Figures\eOPN\ThalamusInactivation\M1_neuralDynamics_Th_inactivation.mat')
neuralDynamics_All_M1_th = neuralDynamics_All;
%%
dynamics = neuralDynamics_All_M1_th;
dimension = 1;
speedTotBaseline = [];
speedToteOPN = [];
rtbaseline = [];
rteopn = [];


for n = 1:length(dynamics)
    % Segment opto and baseline trials
%     [IntanBehaviourBaseline,IntanBehaviourOpto, Waves, WavesOpto] = separateOptoTrials(dynamics(n).IntanBehaviour,dynamics(n).IntanBehaviour.parameters);
%     baselineId = 1:length(IntanBehaviourBaseline.cueHitTrace);
%     eOPNId = length(IntanBehaviourBaseline.cueHitTrace)+1:length(dynamics(n).IntanBehaviour.cueHitTrace);
%     assert(length(eOPNId)==length(IntanBehaviourOpto.cueHitTrace))
    rtbaseline{n} = 0.169;
    rteopn{n} = 0.176;
    speedTotBaseline{n} = squeeze(dynamics{n}.hitbaseline.speed.speed(dimension,2:end,:));
    speedToteOPN{n} = squeeze(dynamics{n}.hitopto.speed.speed(dimension,2:end,:));
    disp(['Calculating session ' num2str(n) '...'])
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
box off,set(gca,'tickdir','out','fontsize',14),axis square,xlim([-500 1500]),ylim([0.01 0.08])
xline(0, '--r', 'Cue');
% xline(mean(rtbaseline)*1000, '--g', 'MI');
% xline(nanmean(rteopn>.400)*1000, '--b', 'MI');
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
significant_pre_cue = mean(p_values(1:50));
disp(['Significant val from pre cue ', num2str((significant_pre_cue))]);
significant_cue_mov = mean(p_values(75:80));
disp(['Significant val from cue to MI ', num2str((significant_cue_mov))]);
significant_mov_post= mean(p_values(85:95));
disp(['Significant val from MI to post ', num2str((significant_mov_post))])

% [p_values, obs_diff, perm_diffs] = paired_signflip_perm(speedTotBaselineSession', speedTotCooledSession', n_permutations);
% significant_cue_mov = mean(p_values);
% disp(['Significant val ', num2str((significant_cue_mov))]);

%%
% % M1 Dynamics for M1 inactivation 
load('\\10.165.57.13\Sutter_backup\Hammad\Ephys\LeverTask\Data_for_Figures\eOPN\M1Inactivation\eOPN_M1_NeuralDynamics_Batch3.mat')
neuralDynamics_All_M1_M1_inact = neuralDynamics_All;
%%
dynamics = neuralDynamics_All_M1_M1_inact;
dimension = 1;
speedTotBaseline = [];
speedToteOPN = [];
rtbaseline = [];
rteopn = [];


for n = 1:length(dynamics)
    if isfield(dynamics{n},'hitbaseline')
        % Segment opto and baseline trials
        %     [IntanBehaviourBaseline,IntanBehaviourOpto, Waves, WavesOpto] = separateOptoTrials(dynamics(n).IntanBehaviour,dynamics(n).IntanBehaviour.parameters);
        %     baselineId = 1:length(IntanBehaviourBaseline.cueHitTrace);
        %     eOPNId = length(IntanBehaviourBaseline.cueHitTrace)+1:length(dynamics(n).IntanBehaviour.cueHitTrace);
        %     assert(length(eOPNId)==length(IntanBehaviourOpto.cueHitTrace))
        rtbaseline{n} = 0.169;
        rteopn{n} = 0.176;
        speedTotBaseline{n} = squeeze(dynamics{n}.hitbaseline.speed.speed(dimension,2:end,:));
%         neuralDynamics_All{n}.hitopto.speed.speed(dimension,2:end,:) = neuralDynamics_All{n}.hitopto.speed.speed(dimension,2:end,:)-0.0066;
        speedToteOPN{n} = squeeze(dynamics{n}.hitopto.speed.speed(dimension,2:end,:));
        disp(['Calculating session ' num2str(n) '...'])
    end
end


speedTotBaseline = horzcat(speedTotBaseline{:});
speedToteOPN = speedToteOPN(~cellfun(@isempty, speedToteOPN));
speedToteOPN = horzcat(speedToteOPN{:});
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
box off,set(gca,'tickdir','out','fontsize',14),axis square,xlim([-500 1500]),ylim([0.01 0.08])
xline(0, '--r', 'Cue');
% xline(mean(rtbaseline)*1000, '--g', 'MI');
% xline(nanmean(rteopn(rteopn>.400))*1000, '--b', 'MI');
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
significant_cue_mov = mean(p_values(77:90));
disp(['Significant val from cue to MI ', num2str((significant_cue_mov))]);
significant_mov_post= mean(p_values(90:110));
disp(['Significant val from MI to post ', num2str((significant_mov_post))]);

% Add values to the current figure/axes
annotationText = sprintf([ ...
    'Mean p-values\n' ...
    'Pre-cue: %.3g\n' ...
    'Cue to MI: %.3g\n' ...
    'MI to post: %.3g'], ...
    significant_pre_cue, significant_cue_mov, significant_mov_post);

text(1500, 0.077, annotationText, ...
    'FontSize', 9, ...
    'VerticalAlignment', 'top', ...
    'BackgroundColor', 'w', ...
    'EdgeColor', [0.6 0.6 0.6], ...
    'Margin', 5);

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