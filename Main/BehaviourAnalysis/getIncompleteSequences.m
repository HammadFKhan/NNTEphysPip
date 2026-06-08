function IntanBehaviour = getIncompleteSequences(IntanBehaviour)
%% Identify incorrect movement sequences from pull and reward timestamps
% Inputs expected in workspace:
%   IntanBehaviour.pullCount   -> 1D array of pull timestamps (ms)
%   rewardIndex               -> 1D array of reward timestamps (ms)

% -----------------------------
% Parameters
% -----------------------------
sequenceBreakMs = 1000;   % New sequence if inter-pull gap > 1000 ms
rewardBufferMs  = 500;    % Reward window extends to 500 ms after last pull

% -----------------------------
% Pull inputs
% -----------------------------
rewardIndex = find(diff(IntanBehaviour.rewardTrace)==1)+1;

pullTimes = IntanBehaviour.pullCount(:);
rewardTimes = rewardIndex(:);

% Remove NaNs if present
pullTimes = pullTimes(~isnan(pullTimes));
rewardTimes = rewardTimes(~isnan(rewardTimes));

% Sort to ensure monotonic order
pullTimes = sort(pullTimes);
rewardTimes = sort(rewardTimes);

% Optional safety check
if isempty(pullTimes)
    error('IntanBehaviour.pullCount is empty. No pull events found.');
end

% -----------------------------
% Sequence segmentation
% -----------------------------
ipi = diff(pullTimes);  % inter-pull intervals
sequenceStartIdx = [1; find(ipi > sequenceBreakMs) + 1];
sequenceEndIdx   = [sequenceStartIdx(2:end) - 1; numel(pullTimes)];
nSequences = numel(sequenceStartIdx);

% -----------------------------
% Preallocate outputs
% -----------------------------
sequenceID       = (1:nSequences).';
startTime        = zeros(nSequences,1);
endTime          = zeros(nSequences,1);
rewardWindowEnd  = zeros(nSequences,1);
nPulls           = zeros(nSequences,1);
label            = strings(nSequences,1);
rewardDetected   = false(nSequences,1);
nRewardsInWindow = zeros(nSequences,1);
isIncorrect      = false(nSequences,1);
isMismatch       = false(nSequences,1);

% Optional detailed bookkeeping
pullIndexStart   = zeros(nSequences,1);
pullIndexEnd     = zeros(nSequences,1);

% -----------------------------
% Sequence-wise analysis
% -----------------------------
for s = 1:nSequences
    idx1 = sequenceStartIdx(s);
    idx2 = sequenceEndIdx(s);

    seqPulls = pullTimes(idx1:idx2);

    pullIndexStart(s) = pullTimes(idx1);
    pullIndexEnd(s)   = pullTimes(idx2);

    startTime(s)       = seqPulls(1);
    endTime(s)         = seqPulls(end);
    rewardWindowEnd(s) = seqPulls(end) + rewardBufferMs;
    nPulls(s)          = numel(seqPulls);

    % Classify sequence by number of pulls
    if nPulls(s) == 1
        label(s) = "single";
    elseif nPulls(s) == 2
        label(s) = "double";
    elseif nPulls(s) == 3
        label(s) = "target";
    else
        label(s) = "extended_bout";
    end

    % Reward check: any reward between first pull and 500 ms after last pull
    rewardMask = rewardTimes >= startTime(s) & rewardTimes <= rewardWindowEnd(s);
    nRewardsInWindow(s) = sum(rewardMask);
    rewardDetected(s) = any(rewardMask);

    % Flagging logic
    % Incorrect: single/double with NO reward
    % Mismatch:  single/double WITH reward
    if (label(s) == "single" || label(s) == "double")
        if rewardDetected(s)
            isMismatch(s) = true;
        else
            isIncorrect(s) = true;
        end
    end
end

% -----------------------------
% Full summary as structure array
% -----------------------------
sequenceSummary = repmat(struct( ...
    'sequenceID', [], ...
    'pullIndexStart', [], ...
    'pullIndexEnd', [], ...
    'startTime', [], ...
    'endTime', [], ...
    'rewardWindowEnd', [], ...
    'nPulls', [], ...
    'label', "", ...
    'rewardDetected', false, ...
    'nRewardsInWindow', [], ...
    'isIncorrect', false, ...
    'isMismatch', false), nSequences, 1);

for s = 1:nSequences
    sequenceSummary(s).sequenceID       = sequenceID(s);
    sequenceSummary(s).pullIndexStart   = pullIndexStart(s);
    sequenceSummary(s).pullIndexEnd     = pullIndexEnd(s);
    sequenceSummary(s).startTime        = startTime(s);
    sequenceSummary(s).endTime          = endTime(s);
    sequenceSummary(s).rewardWindowEnd  = rewardWindowEnd(s);
    sequenceSummary(s).nPulls           = nPulls(s);
    sequenceSummary(s).label            = label(s);
    sequenceSummary(s).rewardDetected   = rewardDetected(s);
    sequenceSummary(s).nRewardsInWindow = nRewardsInWindow(s);
    sequenceSummary(s).isIncorrect      = isIncorrect(s);
    sequenceSummary(s).isMismatch       = isMismatch(s);
end
% -----------------------------
% Final incorrect-sequence structure
% Excludes mismatches by definition
% -----------------------------
incorrectSequenceData = sequenceSummary([sequenceSummary.isIncorrect]);

% -----------------------------
% Optional mismatch structure for QC
% -----------------------------
mismatchSequenceData = sequenceSummary([sequenceSummary.isMismatch]);

% -----------------------------
% Optional compact summary stats
% -----------------------------
allLabels = string({sequenceSummary.label});

summaryStats = struct();
summaryStats.totalSequences     = numel(sequenceSummary);
summaryStats.nSingle            = sum(allLabels == "single");
summaryStats.nDouble            = sum(allLabels == "double");
summaryStats.nTarget            = sum(allLabels == "target");
summaryStats.nExtendedBout      = sum(allLabels == "extended_bout");
summaryStats.nIncorrect         = numel(incorrectSequenceData);
summaryStats.nMismatch          = numel(mismatchSequenceData);
summaryStats.nRewardedSequences = sum([sequenceSummary.rewardDetected]);

% -----------------------------
% Display
% -----------------------------
disp('Full sequence summary:');
disp(sequenceSummary);

disp('Incorrect sequences only:');
disp(incorrectSequenceData);

disp('Mismatch sequences (rewarded single/double, excluded from incorrect data):');
disp(mismatchSequenceData);

disp('Summary stats:');
disp(summaryStats);
%% Now allocate LFP index and time of the data into an incomplete sequence structure
%% Build IntanBehaviour.incompleteSqTrace from incorrectSequenceData
% Assumes incorrectSequenceData is a structure array with fields such as:
%   sequenceID, nPulls, label, firstPullIndex, lastPullIndex
%
% Also assumes:
%   IntanBehaviour.pullCount contains sample indices
%   IntanBehaviour.time is aligned to those indices
%   IntanBehaviour.leverTrace is indexed in the same sample space

parameters = IntanBehaviour.parameters;

nBefore = round(parameters.windowBeforeMI * parameters.Fs);
nAfter  = round(parameters.windowAfterMI  * parameters.Fs);

% Preallocate to max possible size, then trim skipped entries
nCandidates = numel(incorrectSequenceData);

IntanBehaviour.incompleteSqTrace = repmat(struct( ...
    'trace', [], ...
    'time', [], ...
    'LFPIndex', [], ...
    'LFPtime', [], ...
    'pullCount', [], ...
    'sequenceID', [], ...
    'nPulls', [], ...
    'label', "", ...
    'firstPullIndex', [], ...
    'lastPullIndex', []), nCandidates, 1);

keepCount = 0;

for i = 1:nCandidates

    % First pull is the alignment anchor
    firstPullIdx = incorrectSequenceData(i).startTime;
    lastPullIdx  = incorrectSequenceData(i).endTime;

    winStart = firstPullIdx - nBefore;
    winEnd   = firstPullIdx + nAfter;

    % Skip if window exceeds recording bounds
    if winStart < 1 || winEnd > numel(IntanBehaviour.leverTrace) || winEnd > numel(IntanBehaviour.time)
        disp("error on window")
        continue
    end

    keepCount = keepCount + 1;

    thisIdx = (winStart:winEnd).';

    % Extract trace and absolute time/index vectors
    IntanBehaviour.incompleteSqTrace(keepCount).trace    = IntanBehaviour.leverTrace(thisIdx).';
    IntanBehaviour.incompleteSqTrace(keepCount).LFPIndex = thisIdx;
    IntanBehaviour.incompleteSqTrace(keepCount).LFPtime  = IntanBehaviour.time(thisIdx).';

    % Relative time vector, aligned to movementMI window / first pull
    IntanBehaviour.incompleteSqTrace(keepCount).time = ...
        ((0:numel(thisIdx)-1)' ./ parameters.Fs) - parameters.windowBeforeMI;

    % Pulls belonging to this incomplete sequence only, relative to MI window start
    seqPullsAbs = IntanBehaviour.pullCount( ...
        IntanBehaviour.pullCount >= firstPullIdx & ...
        IntanBehaviour.pullCount <= lastPullIdx);

    IntanBehaviour.incompleteSqTrace(keepCount).pullCount = seqPullsAbs - winStart;

    % Metadata copied from incorrectSequenceData
    IntanBehaviour.incompleteSqTrace(keepCount).sequenceID     = incorrectSequenceData(i).sequenceID;
    IntanBehaviour.incompleteSqTrace(keepCount).nPulls         = incorrectSequenceData(i).nPulls;
    IntanBehaviour.incompleteSqTrace(keepCount).label          = incorrectSequenceData(i).label;
    IntanBehaviour.incompleteSqTrace(keepCount).firstPullIndex = firstPullIdx;
    IntanBehaviour.incompleteSqTrace(keepCount).lastPullIndex  = lastPullIdx;
end

% Trim any skipped entries
IntanBehaviour.incompleteSqTrace = IntanBehaviour.incompleteSqTrace(1:keepCount);
