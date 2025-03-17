function [IntanBehaviour] = cleanIntanBehaviour(IntanBehaviour)

%% Rejecting diengage miss trials
[IntanBehaviour] = removeDisengageTrials(IntanBehaviour);

%% Rejecting all trials with RT < 0 and more than 1.5 sec
disp('Rejecting all trials with RT less than 0 and more than 1.5 sec ...');
maxRT = 1.5; % in seconds
 
IntanBehaviour.reactionTime = cell2mat(arrayfun(@(s) s.reactionTime, IntanBehaviour.cueHitTrace, 'UniformOutput', false));
badRT = find(IntanBehaviour.reactionTime<0 | IntanBehaviour.reactionTime>maxRT);
disp(['Rejected  ', string(numel(badRT))]);

%% Rejecting Overlapping Trials
[IntanBehaviour] = removeOverlapTrials(IntanBehaviour);

IntanBehaviour.rejectTrials.Hit = [IntanBehaviour.rejectTrials.Hit badRT];
IntanBehaviour.rejectTrials.Hit = sort(unique(IntanBehaviour.rejectTrials.Hit));

IntanBehaviour.rejectTrials.cueHitTraceReject = IntanBehaviour.cueHitTrace(IntanBehaviour.rejectTrials.Hit);
IntanBehaviour.rejectTrials.hitTraceReject = IntanBehaviour.hitTrace(IntanBehaviour.rejectTrials.Hit);
IntanBehaviour.rejectTrials.MIHitTraceReject = IntanBehaviour.MIHitTrace(IntanBehaviour.rejectTrials.Hit);

IntanBehaviour.rejectTrials.missTraceReject = IntanBehaviour.missTrace(IntanBehaviour.rejectTrials.FA);
IntanBehaviour.rejectTrials.MIFATraceReject = IntanBehaviour.MIFATrace(IntanBehaviour.rejectTrials.FA);

IntanBehaviour.cueHitTrace(IntanBehaviour.rejectTrials.Hit) = [];
IntanBehaviour.hitTrace(IntanBehaviour.rejectTrials.Hit) = [];
IntanBehaviour.MIHitTrace(IntanBehaviour.rejectTrials.Hit) = [];
IntanBehaviour.reactionTime(IntanBehaviour.rejectTrials.Hit) = [];
IntanBehaviour.missTrace(IntanBehaviour.rejectTrials.FA) = [];
IntanBehaviour.MIFATrace(IntanBehaviour.rejectTrials.FA) = [];

IntanBehaviour.cleanFlag = 1;





