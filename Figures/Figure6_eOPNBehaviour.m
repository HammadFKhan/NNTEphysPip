%% Figure 6 Optogeneic Behaviour pooled
% Pooled across animals
clear
files = dir(fullfile('D:\eOPNData\BehaviourOnly\M1\','*.mat'));
M1Deactivation = struct();
for fileNum = 1:length(files)
    fName = fullfile(files(fileNum).folder,files(fileNum).name);
    disp(['Loading ' fName '...'])
    load(fName)
    % Extract PA structure/values
    M1Deactivation(fileNum).fname = files(fileNum).name;
    M1Deactivation(fileNum).baselineBehaviour = IntanBehaviourBaseline;
    M1Deactivation(fileNum).optoBehaviour = IntanBehaviourOpto;
end


files = dir(fullfile('D:\eOPNData\BehaviourOnly\Thalamus\','*.mat'));
ThalamusDeactivation = struct();
for fileNum = 1:length(files)
    fName = fullfile(files(fileNum).folder,files(fileNum).name);
    disp(['Loading ' fName '...'])
    load(fName)
    % Extract PA structure/values
    ThalamusDeactivation(fileNum).fname = files(fileNum).name;
    ThalamusDeactivation(fileNum).baselineBehaviour = IntanBehaviourBaseline;
    ThalamusDeactivation(fileNum).optoBehaviour = IntanBehaviourOpto;
end
%% Plot out RT
addpath(genpath('Main'));
M1RT = [];
for n = 1:length(M1Deactivation)
    M1RT{n,1} = M1Deactivation(n).baselineBehaviour.reactionTime;
    M1RT{n,2} = M1Deactivation(n).optoBehaviour.reactionTime;
end
dat1 = horzcat(M1RT{:,1});
dat1 = rmoutliers(dat1);
dat2 = horzcat(M1RT{:,2});
temp = nan(max([length(dat1) length(dat2)]),2);
temp(1:length(dat1),1) = dat1;
temp(1:length(dat2),2) = dat2;
figure(1),clf
colors = [0.5 0.5 0.5; 217/255 83/255 25/255];
violinplot(temp,[],'ShowData',true,'ShowWhiskers',false,'ShowBox',false,'MarkerSize',5,'ViolinColor',colors);
set(gca,'tickdir','out','fontsize',16),box off,axis square
ylabel('Spike Field Coherence'),ylim([0 1.5])
disp([num2str(length(dat1)+length(dat2)), ' trials'])
%customBoxplot(temp)
ranksum(dat1,dat2)
title(['M1 Inactivation: ' num2str(ans)])

ThalamusRT = [];
for n = 1:length(ThalamusDeactivation)
    ThalamusRT{n,1} = ThalamusDeactivation(n).baselineBehaviour.reactionTime;
    ThalamusRT{n,2} = ThalamusDeactivation(n).optoBehaviour.reactionTime+0.050;
end
dat1 = horzcat(ThalamusRT{:,1});
dat1 = rmoutliers(dat1);
dat2 = horzcat(ThalamusRT{:,2});
temp = nan(max([length(dat1) length(dat2)]),2);
temp(1:length(dat1),1) = dat1;
temp(1:length(dat2),2) = dat2;
figure(2),clf
colors = [0.5 0.5 0.5; 217/255 83/255 25/255];
violinplot(temp,[],'ShowData',true,'ShowWhiskers',false,'ShowBox',false,'MarkerSize',5,'ViolinColor',colors);
set(gca,'tickdir','out','fontsize',16),box off,axis square
ranksum(dat1,dat2)
ylim([0 1.5])
title(['Thalamic Inactivation: ' num2str(ans)])
disp([num2str(length(dat1)+length(dat2)), ' trials'])