function [ripples,i] =  rippleDetection(ripple_signal,timestamps,Fs)
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
    disp('Detection by thresholding failed');
    return
else
    disp(['After detection by thresholding: ' num2str(length(firstPass)) ' events.']);
    
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
    
    if isempty(ripples)
        disp(['No ripples detected on channel ' num2str(size(ripple_signal,2))]);
        i = size(ripple_signal,2)+1;
        return 
    end
    
    %%
    % Optionally, plot results
    figure('Name','Detected Signal')
    plot(timestamps,signal);hold on;
    for j=1:size(ripples,1)
            plot([ripples(j,1) ripples(j,1)],ylim,'g-');
            plot([ripples(j,2) ripples(j,2)],ylim,'k-');
            plot([ripples(j,3) ripples(j,3)],ylim,'r-');
    end
    plot(xlim,[lowThresholdFactor lowThresholdFactor],'k','linestyle','--');
    plot(xlim,[highThresholdFactor highThresholdFactor],'k-');
    axis tight
    box off
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
        print(gcf,'-painters','-depsc', 'Figures/RippleDetect.eps', '-r250');

    
    %% Analyze Ripples
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
    h = hilbert(ripple_signal);
    phase = angle(h);
    amplitude = abs(h);
    unwrapped = unwrap(phase);
    % Compute instantaneous frequency
    frequency = Diff(unwrapped,'smooth',0);
    frequency = frequency/(2*pi);
    
    % Compute ripple map
    [r,i] = Sync([timestamps ripple_signal],ripples(:,2),'durations',durations,'verbose','on');
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
    
    [~,dursort]=sort(data.duration,1,'descend');
    [~,ampsort]=sort(data.peakAmplitude,1,'descend');
    
    x=maps.ripples(dursort,:);
    figure('Name', 'Ripple Map')
    if length(x(:,1)) <= 100
        rippleplot = length(x(:,1));
    else
        rippleplot = 100;
    end
    
    for ii=1:rippleplot
        subplot(10,10,ii)
        plot(x(ii,:))
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        axis off
    end
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
        print(gcf,'-painters','-depsc', 'Figures/Ripple.eps', '-r250');

    figure('Name','SPW-R Data')
    subplot 221
    imagesc(maps.amplitude(ampsort,:))
    title('SPW-R Amplitude: sorted by amplitude')
    subplot 222
    imagesc(maps.amplitude(dursort,:))
    title('SPW-R Amplitude: sorted by duration')
    subplot 223
    imagesc(maps.ripples(ampsort,:))
    title('SPW-R Filtered Signal: sorted by amplitude')
    subplot 224
    imagesc(maps.ripples(dursort,:))
    title('SPW-R Filtered Signal: sorted by duration')
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
        print(gcf,'-painters','-depsc', 'Figures/SPW-R_Amplitude.eps', '-r250');
    
    figure('Name', 'Residuals')
    scatterhist((data.duration*1000),data.peakAmplitude,'kernel','on','Location','SouthWest',...
        'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
        'LineWidth',[2,2,2],'Nbins',[20 100], 'marker','.','markersize',10)
    box off
%     LogScale('x',10)
    xlabel('Duration (ms)'); ylabel('Amplitude (au)')
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
        print(gcf,'-painters','-depsc', 'Figures/Residuals.eps', '-r250');
    
    electrode = 1;
    f = figure;set(f,'name',['Ripple Stats - ' int2str(electrode)]);
    
    subplot(2,2,1);a = gca;hold on;
    plot(((1:nBins)'-ceil(nBins/2))/nBins*diff(durations),maps(electrode).ripples','b');
    
    subplot(2,2,2);
    b = bar(stats(electrode).acg.t,stats(electrode).acg.data);set(b,'FaceColor','r');xlabel('Autocorrelogram');
%      imagesc(stats(electrode).acg.t,stats(electrode).acg.data);colormap(jet);
    %  	b = bar(((0:nCorrBins)-nCorrBins/2)/1000,stats{electrode}.acg.data);xlim([-nCorrBins nCorrBins]/2000);set(b,'FaceColor',[0 0 0]);xlabel('Autocorrelogram');
    
    subplot(2,3,4);a = gca;
    PlotDistribution(data(electrode).peakAmplitude,data(electrode).peakFrequency,'nbins',1000); %,'smooth',5
    axes(a);xlabel(['r=' num2str(stats(electrode).amplitudeFrequency.rho(1,2)) ' p=' num2str(stats(electrode).amplitudeFrequency.p(1,2))]);ylabel('Frequency vs Amplitude');
    
    subplot(2,3,5);a = gca;
    PlotDistribution(data(electrode).duration,data(electrode).peakFrequency,'nbins',1000); %,'smooth',5
    axes(a);xlabel(['r=' num2str(stats(electrode).durationFrequency.rho(1,2)) ' p=' num2str(stats(electrode).durationFrequency.p(1,2))]);ylabel('Frequency vs Duration');
    line([prctile(data(electrode).duration,99.5) prctile(data(electrode).duration,99.5)],[100 200],'color','k')
    
    subplot(2,3,6);a = gca;
    PlotDistribution(data(electrode).duration,data(electrode).peakAmplitude,'nbins',1000); %,'smooth',5
    axes(a);xlabel(['r=' num2str(stats(electrode).durationAmplitude.rho(1,2)) ' p=' num2str(stats(electrode).durationAmplitude.p(1,2))]);ylabel('Amplitude vs Duration');
    
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 3]);...
        print(gcf,'-painters','-depsc', 'Figures/SWR-stats.eps', '-r250');
    
    
end
