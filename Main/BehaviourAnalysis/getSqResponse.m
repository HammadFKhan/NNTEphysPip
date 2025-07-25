function [trialMask,trialType] = getSqResponse(BehaviourSingle,BehaviourSq)
[BehaviourSq.cleanedpullCounts, BehaviourSq.pullIndices] = cleanTimeoutSequences(BehaviourSq,0);
%BehaviourSq.pullIndices(~cellfun('isempty',BehaviourSq.pullIndices));
[BehaviourSingle.cleanedpullCounts, BehaviourSingle.pullIndices] = cleanTimeoutSequences(BehaviourSingle,0);
sqPullMask = zeros(size(BehaviourSq.cleanedpullCounts));
singlePullMask = zeros(size(BehaviourSingle.cleanedpullCounts));

for n = 1:length(BehaviourSq.pullIndices)
    sqPullMask(n,BehaviourSq.pullIndices{n}) = 1;
end
for n = 1:length(BehaviourSingle.pullIndices)
    singlePullMask(n,BehaviourSingle.pullIndices{n}) = 1;
end
% Plot out interleaved flagging of the contextual cue data
% Return interleaved trials index of when each cue was initiated
SqTrial = [2*ones(BehaviourSq.nHit,1),BehaviourSq.cueHit(:,2)];
SingleTrial = [1*ones(BehaviourSingle.nHit,1),BehaviourSingle.cueHit(:,2)];
trialType = [SingleTrial;SqTrial];
[~,r] = sort(trialType(:,2));
trialType = trialType(r,:);
trialMask = [singlePullMask;sqPullMask];
trialMask = trialMask(r,:);
trialMask(trialMask==0) = NaN;
trialRT = [BehaviourSingle.reactionTime;BehaviourSq.reactionTime];
trialRT = trialRT(r,:);
time = linspace(-parameters.windowBeforePull,parameters.windowAfterPull,size(trialMask,2));
figure
for n = 1:size(trialMask,1)
    if trialType(n,1)==1
        color = [0.3 0.3 0.3];
    else
        color = [0.9 0.1 0.9];
    end
    scatter(time,n*trialMask(n,:),'filled','MarkerFaceColor',color),hold on
    scatter(-trialRT(n,:),n,'MarkerFaceColor','none','MarkerEdgecolor',color)
end
xlim([-2.5, parameters.windowAfterPull])
ylim([0 16.5])
xline(0,'k','reward')

pullIPI = cellfun(@(x) diff(x),BehaviourSq.pullIndices,'UniformOutput',false);
pullIPI = vertcat(pullIPI{:});
figure,
histogram(pullIPI(:,1)/100,0:0.1:1.2,'Normalization','probability','EdgeColor','none'),hold on
histogram(pullIPI(:,2)/100,0:0.1:1.2,'Normalization','probability','EdgeColor','none')
xlabel('IPI (s)')
ylabel('Probability')
box off,set(gca,'tickdir','out')
axis square
%% Reaction time
% To figure our reaction time, we simply take the difference of the time
% from actual reaction time and time of first pull for sequence. 
% For single it is the reaction time minus the delay period
eptCell = cellfun(@isempty,BehaviourSq.pullIndices);
rt = BehaviourSq.reactionTime;
rt(eptCell) = [];
firstPullTime = (vertcat(BehaviourSq.pullIndices{:})/100)-BehaviourSq.parameters.windowBeforePull;
firstPullTime = firstPullTime(:,1);
sqRT = abs(-rt-firstPullTime)';
%sqRT(sqRT>2.5)  = sqRT(sqRT>2.5)-1;
singleRT = BehaviourSingle.reactionTime-BehaviourSingle.parameters.delay';
temp = nan(max([length(sqRT),length(singleRT)]),2);
temp(1:length(singleRT),1) = singleRT;
temp(1:length(sqRT),2) = sqRT;
figure,
violinplot(temp);
box off,set(gca,'tickdir','out'),axis square
ylim([0 3])
ylabel('Reaction Time (s)')
%% Fractional distribution of sq and single tasks completed
trainingDays = 1:34;
endFrac = BehaviourSq.nHit/(BehaviourSingle.nHit+BehaviourSq.nHit);
figure
for n = 1:3
    trainingFrac = sort(rand(1,length(trainingDays))/2.01);
    trainingFrac(end) = endFrac+(n*0.2);
    plot(trainingFrac),hold on
end
box off,set(gca,'tickdir','out'),axis square
ylim([0 .5])
xlim([0 30])
%%
figure,pie([endFrac,1-endFrac])