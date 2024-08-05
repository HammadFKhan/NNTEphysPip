function LFP = IntrabetaBurstDetection(signal,Fs)

beta_signal = signal;

timestamps(:,1) = (1:length(beta_signal))/Fs;
LFP.betaBurst = [];
lowThresholdFactor = 1; % Beta envolope must exceed lowThresholdFactor*stdev
highThresholdFactor = 2; % Beta peak must exceed highThresholdFactor*stdev
minInterRippleInterval = 30; % 30ms
minBetaDuration = 40; % 50ms
maxBetaDuration = 200; % 200ms
noise = [];

windowLength = round(11);
data = struct();
stats = struct();
maps = struct();
signal = beta_signal;
squaredSignal = signal.^2;
normalizedSquaredSignal = squaredSignal - mean(squaredSignal)/std(squaredSignal);
% Detect beta periods by thresholding normalized squared signal
thresholded = normalizedSquaredSignal > (lowThresholdFactor*std(squaredSignal));

start = find(diff(thresholded)>0);
stop = find(diff(thresholded)<0);
% Exclude last beta if it is incomplete
if length(stop) == length(start)-1,
    start = start(1:end-1);
end
% Exclude first beta event if it is incomplete
if length(stop)-1 == length(start),
    stop = stop(2:end);
end
% Correct special case when both first and last ripples are incomplete
if start(1) > stop(1),
    stop(1) = [];
    start(end) = [];
end
firstPass = [start,stop];
if isempty(firstPass),
    disp('Detection by thresholding failed')
    return
else
    disp(['After detection by thresholding: ' num2str(length(firstPass)) ' events.']);
end

% Merge beta events if inter-beta period is too short
minInterRippleSamples = minInterRippleInterval/1000*Fs;
secondPass = [];
betaEvent = firstPass(1,:);
for i = 2:size(firstPass,1)
    if firstPass(i,1) - betaEvent(2) < minInterRippleSamples,
        % Merge
        betaEvent = [betaEvent(1) firstPass(i,2)];
    else
        secondPass = [secondPass ; betaEvent];
        betaEvent = firstPass(i,:);
    end
end
secondPass = [secondPass ; betaEvent];
if isempty(secondPass),
    disp('Ripple merge failed');
    return
else
    disp(['After ripple merge: ' num2str(length(secondPass)) ' events.']);
end

% Discard beta events with a peak power < highThresholdFactor
thirdPass = [];
peakNormalizedPower = [];
for i = 1:size(secondPass,1)
    [maxValue,maxIndex] = max(normalizedSquaredSignal([secondPass(i,1):secondPass(i,2)]));
    if maxValue > highThresholdFactor,
        thirdPass = [thirdPass ; secondPass(i,:)];
        peakNormalizedPower = [peakNormalizedPower ; maxValue];
    end
end
if isempty(thirdPass),
    disp('Peak thresholding failed.');
    return
else
    disp(['After peak thresholding: ' num2str(length(thirdPass)) ' events.']);
end

% Detect negative peak position for each beta event
peakPosition = zeros(size(thirdPass,1),1);
for i=1:size(thirdPass,1),
    [minValue,minIndex] = min(signal(thirdPass(i,1):thirdPass(i,2)));
    peakPosition(i) = minIndex + thirdPass(i,1) - 1;
end
% Discard beta events that are way too long
betaBurst = [timestamps(thirdPass(:,1))  timestamps(peakPosition)  timestamps(thirdPass(:,2)) peakNormalizedPower];
duration = betaBurst(:,3)-betaBurst(:,1);
betaBurst(duration>maxBetaDuration/1000,:) = NaN;
disp(['After max duration test: ' num2str(size(betaBurst,1)) ' events.']);

% Discard beta events that are way too short
duration = betaBurst(:,3)-betaBurst(:,1);
betaBurst(duration<minBetaDuration/1000,:) = NaN;
betaBurst = betaBurst((all((~isnan(betaBurst)),2)),:);
disp(['After min duration test: ' num2str(size(betaBurst,1)) ' events.']);

if isempty(betaBurst) || size(betaBurst,1) < 2
    disp(['No events detected on channel']);
    LFP.betaBurst.detectedBeta = [];
end

for i = 1:size(betaBurst,1)
    betaTrace{i} = beta_signal(ceil(Fs*betaBurst(i,1)):ceil((Fs*betaBurst(i,3))));
end
if size(betaBurst,1)==0
    betaTrace = [];
end

LFP.Fs = Fs;
LFP.betaBurst.beta = beta_signal';
LFP.betaBurst.detectedBeta = betaBurst;
LFP.betaBurst.NumDetectedBeta = size(betaBurst,1);
LFP.betaBurst.lowThresholdFactor = lowThresholdFactor;
LFP.betaBurst.highThresholdFactor = highThresholdFactor;
LFP.betaTrace = betaTrace;
