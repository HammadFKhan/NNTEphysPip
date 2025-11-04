function [c,allTrials] = sort_hit_effort(IntanBehaviour)
% sort and concatenate hit and effort trials in order
effortPullTime = nan(1,length(IntanBehaviour.effortperturbTrace));
hitTrialPullTime = nan(1,length(IntanBehaviour.MIHitTrace));
for n = 1:length(IntanBehaviour.hitTrace)
    trueTime = IntanBehaviour.MIHitTrace(n).LFPIndex;
    firstPull = IntanBehaviour.MIHitTrace(n).pullCount(1);
    hitTrialPullTime(n) = trueTime(firstPull);
end
for n = 1:length(IntanBehaviour.effortperturbTrace)
    trueTime = IntanBehaviour.effortperturbTrace(n).LFPIndex;
    firstPull = IntanBehaviour.effortperturbTrace(n).pullCount(1);
    effortPullTime(n) = trueTime(firstPull);
end

% Find the index for shared trials (ie. when a breakout trial was marked
% out)
% It also is required because of some bugs we have....
BO_trials = intersect(hitTrialPullTime,effortPullTime);
if ~isempty(BO_trials)
    for BO = 1:length(BO_trials)
        trialID = find(hitTrialPullTime==BO_trials(BO));
        IntanBehaviour.hitTrace(trialID).effortFlag = 1;
    end
    %assert(length(BO_trials)==length(vertcat(IntanBehaviour.effortperturbTrace.rewardFlag)))
end

allTrials = [hitTrialPullTime,effortPullTime];
allTrials(2,:) = [ones(1,length(hitTrialPullTime)),zeros(1,length(effortPullTime))];
allTrials(3,:) = [1:length(hitTrialPullTime),1:length(effortPullTime)];
[~,c] = sort(allTrials(1,:)); %Sort chronologically

allTrials = allTrials(:,c);
end