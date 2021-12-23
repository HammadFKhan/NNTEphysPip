% Beta analysis
function betaAnalysis(LFP)
detectedBeta = LFP.betaBurst.detectedBeta;
window = LFP.betaBurst.window;
Fs = LFP.downSampleFreq;
if length(detectedBeta)>1
    detectedBeta = vertcat(detectedBeta{:});
    trialFlag = 1;
else
    detectedBeta = cell2mat(detectedBeta);
end

betaDuration = (detectedBeta(:,3)-detectedBeta(:,1)).*1000;
betaAmplitude = detectedBeta(:,4);
CSDoutput = [];
%CSD during velocity window
for i = 1:size(window,1)
    t = LFP.LFP(:,window(i,1)*1024:window(i,2)*1024);
    figure,CSDoutput(:,:,i)  = CSD(t'/1E6,1024,2E-5);
end
mCSD = mean(CSDoutput,3);
figure,imagesc(0:1/Fs:2049,1:64,mCSD'),colormap(jet)
% CSD during beta triggered event
%%
CSDoutput = [];
for i = 1:size(detectedBeta,1)
    t = LFP.LFP(:,detectedBeta(i,1)*Fs:detectedBeta(i,3)*Fs);
    figure,CSDoutput(:,:,i)  = CSD(t'/1E6,1024,2E-5);
end
mCSD = mean(CSDoutput,3);
figure,imagesc(interp2(mCSD',2)),colormap(jet)
%% Plot Beta Events
plt = 1;
if plt ==1
    for i = 1:size(detectedBeta,1)
        eventdiff(i) = detectedBeta(i,3)-detectedBeta(i,1);
    end
    
    figure('Name', 'Residuals')
    scatterhist(betaDuration,betaAmplitude,'kernel','on','Location','SouthWest',...
        'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
        'LineWidth',[2,2,2],'Nbins',[20 100], 'marker','.','markersize',10)
    box off

    % Create array for beta events
    maxdiff = max(eventdiff)*Fs;
    peakWin = ceil(maxdiff/2);
    % betaEvents = zeros(size(detectedBeta{1},1),maxdiff);
    
    % Plot beta events
     for i = 1:size(detectedBeta,1)
            pos1 = detectedBeta(i,1)*Fs;
            pos2 = detectedBeta(i,3)*Fs;
        if trialFlag == 1
            peak = detectedBeta(i,2)*Fs;
            peakAlign(i,:) = LFP.beta_band(peak-peakWin:peak+peakWin);
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
    figure,
    plot(peakAlign','Color',[0.5 0.5 0.5 0.5]); hold on
    plot(mean(peakAlign,1),'k','lineWidth',2)
    xlim([0 maxdiff])
    
    lineError(1:257,peakAlign,'std') % Pass std or ste for error plotting
    
    [wavelet,f] = cwt(mean(peakAlign,1),1024,'FrequencyLimits',[1 30]);
    figure,imagesc((-maxdiff/2):(maxdiff/2),f,abs(wavelet)),colormap(jet)
end
