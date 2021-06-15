function Ripples = rippleAnalysis(filtData,Ripples,Spikes,Fs)
%% Ripple Analysis
if nargin < 4 || strcmp(Fs,'')
    Fs = 20000;
    disp(['Sampling rate set at ' num2str(Fs) ' Hz for ripple analysis']);
end
data = struct();
stats = struct();
maps = struct();
timestamps = filtData.ripple.timestamps;
rippleChannel = 0;
H = waitbar(0,'Analyzing Ripple Waveforms...');
for j = 1:size(Ripples.detectedripples,2)
    if isempty(Ripples.detectedripples{j})== 0
        rippleChannel = rippleChannel+1;
        ripples = Ripples.detectedripples{j};
        signal = filtData.ripple.data(:,j);
        durations = [-0.025 0.025];
        nCorrBins = 50;
        %  corrBinSize = 400;
        corrDuration = 20;
        corrBinSize = 0.01;
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
        
        Ripples(rippleChannel).channelsDetected(rippleChannel) = j;
        Ripples(rippleChannel).ripples = ripples;
        Ripples(rippleChannel).signal = signal;
        Ripples(rippleChannel).timestamps = timestamps;
        %         Ripples(rippleChannel).lowThresholdFactor = lowThresholdFactor;
        %         Ripples(rippleChannel).highThresholdFactor = highThresholdFactor;
        %         Ripples(rippleChannel).nBins = nBins;
        Ripples(rippleChannel).durations = durations;
        Ripples(rippleChannel).stats = stats;
        Ripples(rippleChannel).maps = maps;
        Ripples(rippleChannel).data = data;
        Ripples(rippleChannel).nBins = nBins;
        
        %% SWR onset
        LFPData = filtData.lowpassData';
        rawData = filtData.rawData(j,:)';
        count = 1;
        for ii = 1:size(Ripples.ripples,1)
            start = Fs*(Ripples.ripples(ii,2)-.125);
            stop = Fs*(Ripples.ripples(ii,2)+.125);
            if start > 0 && stop < length(LFPData)
                findDataRaw(:,count) = rawData(start:stop,1);
%                 findSpikes(:,:,count) = Spikes.binary(:,start:stop);
                findDataLFP(:,:,count) = LFPData(start:stop,:);
                [cfs(:,:,count),f] = cwt(findDataRaw(:,count),Fs,'FrequencyLimits',[0.1 600]); 
                count = count+1;
            else
              count = count;  
            end
        end
        Ripples(rippleChannel).rippleOnset.LFP = findDataLFP;
        Ripples(rippleChannel).rippleOnset.SWR = findDataRaw;
%         Ripples(rippleChannel).rippleOnset.PeriStim = findSpikes;
        Ripples(rippleChannel).rippleOnset.cfs = cfs;
        Ripples(rippleChannel).rippleOnset.f = f;
        waitbar(j/size(Ripples.detectedripples,2),H)
        
%         PeriStimt = sum(Ripples.rippleOnset.PeriStim,3);
%         [vectorized,~] = cosine_similarity(PeriStimt(:,1:5000),50);
%         correlation = abs(corr(vectorized));
%         correlation(isnan(correlation)) = 0;
    else
        continue
    end
end

close(H)



