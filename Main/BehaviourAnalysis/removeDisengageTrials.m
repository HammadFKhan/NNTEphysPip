function [IntanBehaviour,Waves] = removeDisengageTrials(IntanBehaviour,Waves)

if exist('Waves','var')
    wavesFlag = 1;
else
    wavesFlag = 0;
end

disp("Rejecting extra miss trials at the end");
extraMisses = 40;
lastHitIndex = IntanBehaviour.cueHitTrace(end).LFPIndex(end);
firstHitIndex = IntanBehaviour.cueHitTrace(1).LFPIndex(1);
IntanBehaviour.rejectTrials.cueMiss = [];
lastMiss = 1;
firstMiss = 0;
for i=1:size(IntanBehaviour.cueMissTrace,2)
    if IntanBehaviour.cueMissTrace(i).LFPIndex(end) < firstHitIndex
        firstMiss = i;
    end
    if IntanBehaviour.cueMissTrace(i).LFPIndex(end) > lastHitIndex
        lastMiss = i;
        break;
    end
end

if firstMiss == 0
    disp('No Miss trials rejected at the start of the trials');
else
    disp(['Rejected first ', string(firstMiss)]);
    IntanBehaviour.rejectTrials.cueMiss = [IntanBehaviour.rejectTrials.cueMiss 1:1:firstMiss];
end

if lastMiss+extraMisses > size(IntanBehaviour.cueMissTrace,2)
    lastMiss = size(IntanBehaviour.cueMissTrace,2);
else
    lastMiss = lastMiss+extraMisses;
    IntanBehaviour.rejectTrials.cueMiss = [IntanBehaviour.rejectTrials.cueMiss lastMiss+1:1:size(IntanBehaviour.cueMissTrace,2)];
end
disp(['Rejected last ', string(size(IntanBehaviour.cueMissTrace,2)-lastMiss)]);



IntanBehaviour.rejectTrials.cueMissTraceRejected = IntanBehaviour.cueMissTrace(IntanBehaviour.rejectTrials.cueMiss);

IntanBehaviour.cueMissTrace(IntanBehaviour.rejectTrials.cueMiss) = [];

if wavesFlag == 1
    Waves.wavesMiss = Waves.wavesMiss(firstMiss+1:lastMiss-1);
end

end

