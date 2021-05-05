function [rippleData,data,stats,maps] =  rippleDetection(ripple_signal,timestamps,Fs)
rippleData = struct();
data = struct();
stats = struct();
maps = struct();
keep = [];
lowThresholdFactor = 2; % Ripple envoloppe must exceed lowThresholdFactor*stdev
highThresholdFactor = 5; % Ripple peak must exceed highThresholdFactor*stdev
minInterRippleInterval = 30; % 30ms
minRippleDuration = 20; % 20ms
maxRippleDuration = 100; % 100ms
noise = [];
windowLength = round(11);
signal = ripple_signal;
squaredSignal = signal.^2;
window = ones(windowLength,1)/windowLength;

normalizedSquaredSignal = squaredSignal - mean(squaredSignal)/std(squaredSignal);
% Detect ripple periods by thresholding normalized squared signal
thresholded = normalizedSquaredSignal > lowThresholdFactor;

start = find(diff(thresholded)>0);
stop = find(diff(thresholded)<0);
% Exclude last ripple if it is incomplete
if length(stop) == length(start)-1,
    start = start(1:end-1);
end
% Exclude first ripple if it is incomplete
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

% Merge ripples if inter-ripple period is too short
minInterRippleSamples = minInterRippleInterval/1000*Fs;
secondPass = [];
ripple = firstPass(1,:);
for i = 2:size(firstPass,1)
    if firstPass(i,1) - ripple(2) < minInterRippleSamples,
        % Merge
        ripple = [ripple(1) firstPass(i,2)];
    else
        secondPass = [secondPass ; ripple];
        ripple = firstPass(i,:);
    end
end
secondPass = [secondPass ; ripple];
if isempty(secondPass),
    disp('Ripple merge failed');
    return
else
    disp(['After ripple merge: ' num2str(length(secondPass)) ' events.']);
end

% Discard ripples with a peak power < highThresholdFactor
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

% Detect negative peak position for each ripple
peakPosition = zeros(size(thirdPass,1),1);
for i=1:size(thirdPass,1),
    [minValue,minIndex] = min(signal(thirdPass(i,1):thirdPass(i,2)));
    peakPosition(i) = minIndex + thirdPass(i,1) - 1;
end
% Discard ripples that are way too long
ripples = [timestamps(thirdPass(:,1))  timestamps(peakPosition)  timestamps(thirdPass(:,2)) peakNormalizedPower];
duration = ripples(:,3)-ripples(:,1);
ripples(duration>maxRippleDuration/1000,:) = NaN;
disp(['After max duration test: ' num2str(size(ripples,1)) ' events.']);

% Discard ripples that are way too short
duration = ripples(:,3)-ripples(:,1);
ripples(duration<minRippleDuration/1000,:) = NaN;
ripples = ripples((all((~isnan(ripples)),2)),:);
disp(['After min duration test: ' num2str(size(ripples,1)) ' events.']);

if isempty(ripples) || size(ripples,1) < 2
    disp(['No ripples detected on channel ' num2str(size(ripple_signal,2))]);
    return
end
%% Ripple Analysis
durations = [-0.025 0.025];
nCorrBins = 50;
%  corrBinSize = 400;
corrDuration = 20;
corrBinSize = 0.01;
signal = ripple_signal; 
nBins = floor(Fs*diff(durations)/2)*2+1; % must be odd
nHalfCenterBins = 3;
centerBin = ceil(nBins/2);
centerBins = centerBin-nHalfCenterBins:centerBin+nHalfCenterBins;
idx = ceil((ripples(:,2)-ripples(:,1))*Fs);
% Compute instantaneous phase and amplitude
h = hilbert(signal);
phase = angle(h);
amplitude = abs(h);
unwrapped = unwrap(phase);
% Compute instantaneous frequency
frequency = Diff(unwrapped,'smooth',0);
frequency = frequency/(2*pi);

% Compute ripple map
[r,i] = Sync([timestamps signal],ripples(:,2),'durations',durations,'verbose','on');
maps.ripples = SyncMap(r,i,'durations',durations,'nbins',nBins,'smooth',0);

% Compute frequency Map
[f,i] = Sync([timestamps frequency],ripples(:,2),'durations',durations);
maps.frequency = SyncMap(f,i,'durations',durations,'nbins',nBins,'smooth',0);

% Compute phase map
[p,i] = Sync([timestamps phase],ripples(:,2),'durations',durations);
maps.phase = SyncMap(p,i,'durations',durations,'nbins',nBins,'smooth',0);

% Compute amplitude map
[a,i] = Sync([timestamps amplitude],ripples(:,2),'durations',durations);
maps.amplitude = SyncMap(a,i,'durations',durations,'nbins',nBins,'smooth',0);

idx(idx>length(maps.frequency(1,:))) = length(maps.frequency(1,:));

% Ripple frequency and amplitude at peak
data.peakFrequency = maps.frequency(:,centerBin);
data.peakAmplitude = maps.amplitude(:,centerBin);

% Ripple durations
data.duration = abs(diff(ripples(:,[1 3])'))';

% Autocorrelogram and correlations
%  if nargin > 2,
%  	[stats.acg.data,stats.acg.t] = CCG(ripples(:,2),1,'binSize',corrBinSize,'halfBins',nCorrBins/2);
[stats.acg.data,stats.acg.t] = CCG(ripples(:,2),ones(length(ripples(:,2)),1),'binSize',corrBinSize);
[stats.amplitudeFrequency.rho,stats.amplitudeFrequency.p] = corrcoef(data.peakAmplitude,data.peakFrequency);
[stats.durationFrequency.rho,stats.durationFrequency.p] = corrcoef(data.duration,data.peakFrequency);
[stats.durationAmplitude.rho,stats.durationAmplitude.p] = corrcoef(data.duration,data.peakAmplitude);

rippleData.signal = signal;
rippleData.ripples = ripples;
rippleData.timestamps = timestamps;
rippleData.lowThresholdFactor = lowThresholdFactor;
rippleData.highThresholdFactor = highThresholdFactor;
rippleData.nBins = nBins;
rippleData.durations = durations;
end
