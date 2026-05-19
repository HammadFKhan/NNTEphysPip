%% Setup motion response
% From the motion energy structure we can extract trial dependant
% kinematics
% Segement data as a function of hit trials
% Grab time length of the camera; 

fsCam  = 30;                     % camera sampling rate (Hz)
nFrames = size(motSVD,1);
timeCamera = (0:nFrames-1)' / fsCam;

parts = {'Tongue','Limb','Whiskers','Body','Pupil'};
nTrials = length(IntanBehaviour.hitTrace);

motionData.hitTrial = struct([]);
trialLength = [];               % will store the expected length

for nTrial = 1:nTrials
    trialStart = IntanBehaviour.hitTrace(nTrial).LFPtime(1);
    trialEnd   = IntanBehaviour.hitTrace(nTrial).LFPtime(end);

    % Find nearest camera indices
    [~, startIdx] = min(abs(timeCamera - trialStart));
    [~, stopIdx]  = min(abs(timeCamera - trialEnd));

    if stopIdx < startIdx
        error('Trial %d has stopIdx < startIdx', nTrial);
    end

    idx = startIdx:stopIdx;
    thisLen = numel(idx);

    % --- length consistency check ---
    if isempty(trialLength)
        trialLength = thisLen;   % first trial sets expected lengthdat
    elseif thisLen ~= trialLength
        error('Trial %d has length %d, expected %d', ...
              nTrial, thisLen, trialLength);
    end
    % --------------------------------

    motionData.hitTrial(nTrial).timeCamera = timeCamera(idx);
    motionData.hitTrial(nTrial).trialStart = trialStart;
    motionData.hitTrial(nTrial).trialEnd   = trialEnd;
    motionData.hitTrial(nTrial).startIdx   = startIdx;
    motionData.hitTrial(nTrial).stopIdx    = stopIdx;

    for p = 1:numel(parts)
        fieldName = [parts{p} '_motion'];
        motionData.hitTrial(nTrial).(fieldName) = bodyPartData.(fieldName)(idx);
    end
end

fprintf('All %d trials have equal length: %d frames.\n', nTrials, trialLength);