%% Figure 5 Trajectory wave coupling during cooling
% Pooled across animals
clear
files = dir(fullfile('D:\TrajectoryWavesCouplingCooling\','*.mat'));
CoupledCooling = struct();
for fileNum = 1:length(files)
    fName = fullfile(files(fileNum).folder,files(fileNum).name);
    disp(['Loading ' fName '...'])
    load(fName)
    % Extract PA structure/values
    CoupledCooling(fileNum).fname = files(fileNum).name;
    CoupledCooling(fileNum).baselineDat = baselineDat;
    CoupledCooling(fileNum).coolingDat = coolingDat;
end

%% Plot out the data and do stats

figure(1),clf
close all
% Check across subset of trials
PGDDiff = [];
SpeedDiff = [];
count = 1;
for n = 1:length(CoupledCooling)%[2 3 5]
    for nn = 1:2
        trialIdBase = randperm(size(CoupledCooling(n).baselineDat.PGDTrajectorySpeed,1));
        trialIdBase = trialIdBase(1:floor(length(trialIdBase)*.75));
        trialIdCool = randperm(size(CoupledCooling(n).coolingDat.PGDTrajectorySpeed,1));
        trialIdCool = trialIdCool(1:floor(length(trialIdCool)*.75));
        PGDDiff(count,:) = abs((mean(CoupledCooling(n).baselineDat.PGDTrajectorySpeed(trialIdBase,:))-mean(CoupledCooling(n).coolingDat.PGDTrajectorySpeed(trialIdCool,:))))/2;
        SpeedDiff(count,:) = abs((mean(CoupledCooling(n).baselineDat.SpeedTrajectorySpeed(trialIdBase,:))-mean(CoupledCooling(n).coolingDat.SpeedTrajectorySpeed(trialIdCool,:))))/2;
        count = count+1;
    end
end
%%
load('D:\TrajectoryWavesCouplingCooling\pooledData.mat')

close all
figure,
subplot(121),errorbar(1:3, mean(PGDDiff,1),std(PGDDiff)/sqrt(5))
xlim([0.5 3.5]),ylim([0 0.4])

axis square,box off,set(gca,'tickdir','out','fontsize',14),axis square,ylabel('PGD Coupling difference')
subplot(122),errorbar(1:3, mean(SpeedDiff,1),std(SpeedDiff)/sqrt(5)),
xlim([0.5 3.5])
axis square,box off,set(gca,'tickdir','out','fontsize',14),axis square,ylabel('Speed Coupling difference')
ylim([0 0.4])

[~,~,stats1] = anova1(PGDDiff);
[~,~,stats2] = anova1(SpeedDiff);
c1 = multcompare(stats1);
c2 = multcompare(stats2);