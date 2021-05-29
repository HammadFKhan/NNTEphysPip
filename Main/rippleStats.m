function [rippleData,data,stats,maps] = rippleStats(ripples,ripple_signal,timestamps,Fs)
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
