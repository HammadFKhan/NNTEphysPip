function betaGroup = groupBetaBurstDetection(LFP,betaTrials,window,Fs)
if nargin<2
    beta_signal = LFP.beta_band';
    timestamps(:,1) = (1:length(beta_signal))/Fs; % Timepoint for each sample
    window = [];
    LFP.betaBurst = [];
else
    beta_signal = squeeze(betaTrials);
    trialFlag = 1;
    betaBurst = [];
end

lowThresholdFactor = 1; % Beta envolope must exceed lowThresholdFactor*stdev
highThresholdFactor = 1; % Beta peak must exceed highThresholdFactor*stdev
minInterRippleInterval = 30; % 30ms
minBetaDuration = 40; % 50ms
maxBetaDuration = 250; % 200ms
noise = [];

windowLength = round(11);
betaGroup.betaBurst.lowThresholdFactor = lowThresholdFactor;
betaGroup.betaBurst.highThresholdFactor = highThresholdFactor;
betaGroup.betaBurst.window = window;
betaGroup.downSampleFreq = Fs;
betaGroup.beta_band = LFP;

H = waitbar(0,'Detecting Beta Events...');
for idx = 1:size(beta_signal,2)
    betaGroup.betaBurst.detectedBeta{idx} = [];
    betaGroup.betaBurst.NumDetectedBeta(idx,1) =0;
    waitbar(idx/size(beta_signal,1),H)
    data = struct();
    stats = struct();
    maps = struct();
    if trialFlag
        timestamps = window(idx,1):1/Fs:window(idx,2);
        timestamps = timestamps';
    end
    normByTrial = 1;
    % super unoptimized method in creating the thresholding as a function
    % of all trials instead of each trial (TODO: check differences but I
    % think not norming by trials is better....
    if normByTrial
        disp('Thresholding by trial...')
        signal = beta_signal(:,idx);
        squaredSignal = signal.^2;
        normalizedSquaredSignal = (squaredSignal - mean(squaredSignal))/std(squaredSignal);
        % Detect beta periods by thresholding normalized squared signal
        thresholded = normalizedSquaredSignal > lowThresholdFactor;
    else 
        disp('Thresholding by entire behavior state...')
        signal = beta_signal(:); %collapes all trials and then threshold instead
        squaredSignal = signal.^2;
        normalizedSquaredSignal = (squaredSignal - mean(squaredSignal))/std(squaredSignal);
        % Detect beta periods by thresholding normalized squared signal
        thresholdedAll = normalizedSquaredSignal > lowThresholdFactor;
        signal = beta_signal(:,idx); %rebuild the signal
        thresholdedAll = reshape(thresholdedAll,size(beta_signal,1),size(beta_signal,2)); %rebuild the threshold 
        thresholded = thresholdedAll(:,idx);                                             % and then take the correct index in the for loop
    end
    
    start = find(diff(thresholded)>0);
    stop = find(diff(thresholded)<0);
    % Check to see if there are any signals thersholded
    if isempty(start) || isempty(stop)
        disp('No candidates detected')
        
    else
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
        peakPower = [];
        for i = 1:size(secondPass,1)
            [maxValue,maxIndex] = max(normalizedSquaredSignal([secondPass(i,1):secondPass(i,2)]));
            [maxValue1,maxIndex2] = max(squaredSignal([secondPass(i,1):secondPass(i,2)]));
            if maxValue > highThresholdFactor,
                thirdPass = [thirdPass ; secondPass(i,:)];
                peakNormalizedPower = [peakNormalizedPower ; maxValue];
                peakPower = [peakPower ; maxValue1];
            end
        end
        if isempty(thirdPass),
            disp('Peak thresholding failed.');
            betaGroup.betaBurst.detectedBeta{idx} = [];
        else
            disp(['After peak thresholding: ' num2str(length(thirdPass)) ' events.'])
            
            % Detect negative peak position for each beta event
            peakPosition = zeros(size(thirdPass,1),1);
            for i=1:size(thirdPass,1),
                [minValue,minIndex] = min(signal(thirdPass(i,1):thirdPass(i,2)));
                peakPosition(i) = minIndex + thirdPass(i,1) - 1;
            end
            % Discard beta events that are way too long
            betaBurst = [timestamps(thirdPass(:,1))  timestamps(peakPosition)  timestamps(thirdPass(:,2)) peakNormalizedPower peakPower];
            duration = betaBurst(:,3)-betaBurst(:,1);
            betaBurst(duration>maxBetaDuration/1000,:) = NaN;
            disp(['After max duration test: ' num2str(size(betaBurst,1)) ' events.']);
            
            % Discard beta events that are way too short
            duration = betaBurst(:,3)-betaBurst(:,1);
            betaBurst(duration<minBetaDuration/1000,:) = NaN;
            betaBurst = betaBurst((all((~isnan(betaBurst)),2)),:);
            disp(['After min duration test: ' num2str(size(betaBurst,1)) ' events.']);
            
            if isempty(betaBurst) || size(betaBurst,1) < 2
                disp(['No events detected on channel ' num2str(idx)]);
                betaGroup.betaBurst.detectedBeta{idx} = [];
                continue
            end
            
            %     betaGroup.LFP = LFP.LFP;
            betaGroup.betaBurst.detectedBeta{idx} = betaBurst;
            betaGroup.betaBurst.NumDetectedBeta(idx,1) = size(betaBurst,1);
        end
    end
end
close(H)
