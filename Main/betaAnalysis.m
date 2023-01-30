% Beta analysis
function [peakAlign,csdPeakAlign,norm,stats] = betaAnalysis(LFP,allLFP)
if nargin<2
    allLFP = [];
    csd = [];
end
detectedBeta = LFP.betaBurst.detectedBeta;
% window = LFP.betaBurst.window;
Fs = LFP.downSampleFreq;

if length(detectedBeta)>1
    detectedBeta = vertcat(detectedBeta{:});
    trialFlag = 1;
else
    detectedBeta = cell2mat(detectedBeta);
end
% Check data validity
if sum(LFP.betaBurst.NumDetectedBeta)==0
    peakAlign = [];
    norm = [];
    csdPeakAlign = [];
    stats.betaDuration = [];
    stats.betaAmplitudeNorm =[];
    stats.betaAmplitudeErrorNorm = [];
    stats.betaAmplitude = [];
    stats.betaAmplitudeError = [];
    stats.betaAmplitudeNum = [];
    stats.betaER = [];
    stats.betaAmplitudeSingle = [];
    stats.betaAmplitudeRaw = [];


    %FanoFactor calculation
    stats.DurationFanoF = [];
    stats.DurationCVsquared = [];
    stats.AmplitudeFanoF = [];
    stats.AmplitudeCVsquared = [];
    stats.freqPks = [];
    disp('No statistics to analyze')
else
    
    stats.betaDuration = (detectedBeta(:,3)-detectedBeta(:,1)).*1000;
    stats.betaAmplitudeNorm = mean(detectedBeta(:,4));
    stats.betaAmplitudeErrorNorm = std(detectedBeta(:,4));
    stats.betaAmplitude = mean(detectedBeta(:,5));
    stats.betaAmplitudeError = std(detectedBeta(:,5));
    stats.betaAmplitudeNum = length(detectedBeta);
    stats.betaER = LFP.betaBurst.NumDetectedBeta(LFP.betaBurst.NumDetectedBeta~=0);
    stats.betaAmplitudeSingle = detectedBeta(:,4);
    stats.betaAmplitudeRaw = detectedBeta(:,5);
    duration = (detectedBeta(:,3)-detectedBeta(:,1)).*1000;
    amplitude = detectedBeta(:,4);
    
    %FanoFactor calculation
    stats.DurationFanoF = std(duration)^2/mean(duration);
    stats.DurationCVsquared = (std(duration)/mean(duration))^2;
    stats.AmplitudeFanoF = std(amplitude)^2/mean(amplitude);
    stats.AmplitudeCVsquared = (std(amplitude)/mean(amplitude))^2;
    
    %% Plot Beta Events
    plt = 1;
    if plt ==1
        for i = 1:size(detectedBeta,1)
            eventdiff(i) = detectedBeta(i,3)-detectedBeta(i,1);
        end
        
        %     figure('Name', 'Residuals')
        %     scatterhist(duration,amplitude,'kernel','on','Location','SouthWest',...
        %         'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
        %         'LineWidth',[2,2,2],'Nbins',[20 100], 'marker','.','markersize',10)
        %     box off
    end
    % Create array for beta events
    maxdiff = 257;
    peakWin = ceil(maxdiff/2);
    % betaEvents = zeros(size(detectedBeta{1},1),maxdiff);
    
    % Plot beta events
    for i = 1:size(detectedBeta,1)
        if trialFlag == 1
            peak = detectedBeta(i,2)*Fs;
            peakAlign(i,:) = LFP.beta_band(peak-peakWin:peak+peakWin);
            csdPeakAlign(:,:,i) = allLFP(:,peak-peakWin:peak+peakWin); %reference to all electrodes
            %             CSDoutput(:,:,i) = zscore(CSD(csdPeakAlign'/1E6,1024,2E-5));
        else
            peakAlign(i,:) = LFP.beta_band(peak-peakWin:peak+peakWin);
        end
        %     originalSignal = LFP.beta_band(1,pos1:pos2);
        %     padding = maxdiff-size(originalSignal,2);
        %     if padding < 0
        %         betaEvents(i,:) = originalSignal(1,1:end+padding);
        %     else
        %         betaEvents(i,:) = horzcat(originalSignal,zeros(1,padding));
        %     end
    end
    %      csd = CSDoutput;
    %     figure,
    %     plot(peakAlign','Color',[0.5 0.5 0.5 0.5]); hold on
    %     plot(mean(peakAlign,1),'k','lineWidth',2)
    %     xlim([0 maxdiff])
    %     figure,lineError(1:size(peakAlign,2),peakAlign,'std') % Pass std or ste for error plotting
    
    %% Beta event profile and frequency event span
    
    [wavelet,f] = cwt(mean(peakAlign,1),1024,'FrequencyLimits',[1 30]);
    blah = abs(wavelet);
    norm = (blah-min(blah,[],'all'))/(max(blah,[],'all')-min(blah,[],'all'));
    %     figure,imagesc((-maxdiff/2):(maxdiff/2),f,norm),colormap(jet),colorbar
    % Calculate the FWHM of each beta event to define freq span
    
    for i = 1:size(peakAlign,1)
        [betaWavelet(:,:,i),betaF] = cwt(peakAlign(i,:), 1024, 'FrequencyLimits',[1 30]);
    end
    %%
    for i = 1:size(betaWavelet,3)
        t = mean(betaWavelet(:,:,i),2);
        t1 = interp1(1:length(t),t,1:0.05:length(t));
        t2(:,i) = smoothdata(abs(t1),'gaussian',50);
        F = interp1(1:length(betaF),betaF,1:0.05:length(betaF));
    end
    
    for i = 1:size(t2,2)
        pks = find(t2(:,i)==max(t2(:,i)));
        fF = sort(F,'ascend'); % need to sort freuquency as ASCENDING!!
        freqPks(i) = fF(pks);
    end
    
    stats.freqPks = freqPks;
end
