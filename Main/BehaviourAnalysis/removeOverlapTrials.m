function [IntanBehaviour] = removeOverlapTrials(IntanBehaviour)

%% Rejects all trials that have overlaps. 
   % 1 - Hits and FA overlap 
   % 2 - Two or more FA overlap

overlapThreshold = 0.5 * IntanBehaviour.parameters.Fs; %time
overlapThresholdUpper = 2.5 * IntanBehaviour.parameters.Fs;

goodTrials.Hit = 1:1:length(IntanBehaviour.MIHitTrace);
goodTrials.HitFA = 1:1:length(IntanBehaviour.MIHitTrace);
goodTrials.FA = 1:1:length(IntanBehaviour.MIFATrace);
goodTrials.FAHit = 1:1:length(IntanBehaviour.MIFATrace);

% timeHit = zeros(size(IntanBehaviour.time));
% timeFA = zeros(size(IntanBehaviour.time));
timeMIHit = zeros(size(IntanBehaviour.time));
timeMIFA = zeros(size(IntanBehaviour.time));

for i=1:size(IntanBehaviour.cueHitTrace,2)
    timeMIHit(IntanBehaviour.MIHitTrace(i).LFPIndex(1):IntanBehaviour.MIHitTrace(i).LFPIndex(end)) = 2;
end

for i=1:size(IntanBehaviour.MIFATrace,2)
    timeMIFA(IntanBehaviour.MIFATrace(i).LFPIndex(1):IntanBehaviour.MIFATrace(i).LFPIndex(end)) = 1;
end

%% Checking for Hit trials
for i=1:size(IntanBehaviour.MIHitTrace,2)
    % Overlap - MIHit and MIFA
    a = zeros(size(IntanBehaviour.time));
    a(IntanBehaviour.MIHitTrace(i).LFPIndex(1):IntanBehaviour.MIHitTrace(i).LFPIndex(end)) = 1;
    if sum(a.*timeMIFA) > overlapThreshold
        if sum(a.*timeMIFA) > overlapThresholdUpper
            goodTrials.HitFA(i) = 0;
            for j=1:size(IntanBehaviour.MIFATrace,2)
                hitTime = IntanBehaviour.cueHitTrace(i).LFPIndex(IntanBehaviour.parameters.windowBeforePull*IntanBehaviour.parameters.Fs);
                if hitTime > IntanBehaviour.MIFATrace(j).LFPIndex(1) && hitTime < IntanBehaviour.MIFATrace(j).LFPIndex(end)
                    goodTrials.FAHit(j) = 0;
                    break;
                end
            end
        else
            goodTrials.Hit(i) = 0;
        end
    end
end

%% Checking for FA trials
currentFA = cell2mat(arrayfun(@(s) s.LFPIndex(end), IntanBehaviour.MIFATrace, 'UniformOutput', false));
nextFA = cell2mat(arrayfun(@(s) s.LFPIndex(1), IntanBehaviour.MIFATrace, 'UniformOutput', false));
currentFA = currentFA(1:end-1);
nextFA = nextFA(2:end);
overlapFA = nextFA-currentFA;

a = find(overlapFA<-overlapThreshold);
a = [a a+1];
a = sort(unique(a));
goodTrials.FA(a) = 0;

%% Rejecting overlapping trials

IntanBehaviour.goodTrials.Hit = nonzeros(goodTrials.Hit);
IntanBehaviour.goodTrials.FA = find(nonzeros(goodTrials.FA));

IntanBehaviour.rejectTrials.Hit = find(goodTrials.Hit == 0);
IntanBehaviour.rejectTrials.FA = find((goodTrials.FA .* goodTrials.FAHit)== 0);
IntanBehaviour.rejectTrials.FAHit = find(goodTrials.FAHit == 0);
IntanBehaviour.rejectTrials.HitFA = find(goodTrials.HitFA == 0);

disp(['Rejected Hits - ', string(numel(goodTrials.Hit)-nnz(goodTrials.Hit))]);
disp(['Rejected FAs - ', string(numel(goodTrials.FA)-nnz(goodTrials.FA.*goodTrials.FAHit))]);

