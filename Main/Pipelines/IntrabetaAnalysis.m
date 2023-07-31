% Beta analysis
function [peakAlign,norm,f,stats] = IntrabetaAnalysis(LFP)

detectedBeta = LFP.betaBurst.detectedBeta;
% window = LFP.betaBurst.window;
Fs = LFP.Fs;

% Check data validity
if sum(LFP.betaBurst.NumDetectedBeta)==0
    peakAlign = [];
    csd = [];
    norm =[];
    f = [];
    stats.betaDuration = 0;
    stats.betaAmplitude = 0;
    stats.betaER = 0;
    disp('No statistics to analyze')
    return
end
stats.betaDuration = (detectedBeta(:,3)-detectedBeta(:,1)).*1000;
stats.betaAmplitudeNorm = mean(detectedBeta(:,4));
stats.betaAmplitudeErrorNorm = std(detectedBeta(:,4));
% stats.betaAmplitude = mean(detectedBeta(:,5));
% stats.betaAmplitudeError = std(detectedBeta(:,5));
stats.betaAmplitudeNum = length(detectedBeta);
stats.betaER = LFP.betaBurst.NumDetectedBeta(LFP.betaBurst.NumDetectedBeta~=0);  
duration = (detectedBeta(:,3)-detectedBeta(:,1)).*1000;
amplitude = detectedBeta(:,4);

%FanoFactor calculation
stats.DurationFanoF = std(duration)^2/mean(duration);
stats.DurationCVsquared = (std(duration)/mean(duration))^2;
stats.AmplitudeFanoF = std(amplitude)^2/mean(amplitude);
stats.AmplitudeCVsquared = (std(amplitude)/mean(amplitude))^2;

CSDoutput = [];
%CSD during velocity window
% for i = 1:size(window,1)
%     t = LFP.LFP(:,window(i,1)*1024:window(i,2)*1024);
%     CSDoutput  = CSD(csdPeakAlign'/1E6,1024,2E-5);
% end
% mCSD = mean(CSDoutput,3);
% figure,imagesc(0:1/Fs:2049,1:64,mCSD'),colormap(jet)
% CSD during beta triggered event
%%
% CSDoutput = [];
% for i = 1:size(detectedBeta,1)
%     t = LFP.LFP(:,detectedBeta(i,1)*Fs:detectedBeta(i,3)*Fs);
%     figure,CSDoutput(:,:,i)  = CSD(t'/1E6,1024,2E-5);
% end
% mCSD = mean(CSDoutput,3);
% figure,imagesc(interp2(mCSD',2)),colormap(jet)
%% Plot Beta Events
plt = 1;
if plt ==1
    for i = 1:size(detectedBeta,1)
        eventdiff(i) = detectedBeta(i,3)-detectedBeta(i,1);
    end
    
    figure('Name', 'Residuals')
    scatterhist(duration,amplitude,'kernel','on','Location','SouthWest',...
        'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
        'LineWidth',[2,2,2],'Nbins',[20 100], 'marker','.','markersize',10)
    box off

    % Create array for beta events
    maxdiff = Fs*0.254;
    peakWin = ceil(maxdiff/2);
    % betaEvents = zeros(size(detectedBeta{1},1),maxdiff);
    
    % Plot beta events
    for i = 1:(size(detectedBeta,1)-1)
                peak = detectedBeta(i,2)*Fs;
                peakAlign(i,:) = LFP.betaBurst.beta(peak-peakWin:peak+peakWin);
%         pos1 = detectedBeta(i,1)*Fs;
%         pos2 = detectedBeta(i,3)*Fs;
%         originalSignal = LFP.betaBurst.beta(1,pos1:pos2);
%         padding = maxdiff-size(originalSignal,2);
%         if padding < 0
%             betaEvents(i,:) = originalSignal(1,1:end+padding);
%         elseif padding>0 && rem(padding,2)==1
%             betaEvents(i,:) = horzcat(zeros(1,floor(padding/2)),originalSignal,zeros(1,ceil(padding/2)));
%         else
%             betaEvents(i,:) = padarray(originalSignal,[0,floor(padding/2)],'both');
%             
%         end
    end
%      csd = CSDoutput;
    figure,
    plot(peakAlign','Color',[0.5 0.5 0.5 0.5]); hold on
    plot(mean(peakAlign,1),'k','lineWidth',2)
    xlim([0 maxdiff])
    
%Beta event profile
    
%     figure,lineError(1:size(peakAlign,2),peakAlign,'std') % Pass std or ste for error plotting
    [wavelet,f] = cwt(smoothdata(mean(peakAlign,1)),Fs,'FrequencyLimits',[1 30]);
    blah = abs(wavelet);
    norm = (blah-min(blah,[],'all'))/(max(blah,[],'all')-min(blah,[],'all'));
    figure,imagesc((-maxdiff/2):(maxdiff/2),f,norm),colormap(jet),colorbar
end
