%% Code for SI beta event detection
% Wavelet of trial # and subsequency filtered Beta LFP
for i = 1
signal = TimeFreq.tfRun.betaLFP(:,1,i);
[wavelet,f] = cwt(signal,1024,'FrequencyLimits',[10 30]);
blah = abs(wavelet);
norm = (blah-min(blah,[],'all'))/(max(blah,[],'all')-min(blah,[],'all'));
figure,imagesc(-1000:1000,f,norm),colormap(jet),colorbar,axis xy,caxis([0.35 .93])
figure,plot(smoothdata(signal,'gaussian',25)),box off
refEvent = betaGroup(1).electrode.betaBurst.detectedBeta{i}(:,2);
refEvent = refEvent-betaGroup(64).electrode.betaBurst.window(i,1);
refEvent = refEvent.*1024;
for ii = 1:length(refEvent)
    xline(refEvent(ii),'r--')
end
end
%% Plot beta amplitude vs duration
detectedBeta = [];
for i = 60:64
detectedBetas = vertcat(betaGroup(i).electrode.betaBurst.detectedBeta{:});
detectedBeta = vertcat(detectedBeta,detectedBetas);
end
figure('Name', 'Residuals')
duration = ((detectedBeta(:,3)-detectedBeta(:,1)).*1000)+10;
amplitude = detectedBeta(:,4);
scatterhist(duration,amplitude,'kernel','on','Location','SouthWest',...
    'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
    'LineWidth',[2,2,2],'Nbins',[20 100], 'marker','.','markersize',10)
box off,hold on
disp(['Beta N: ' num2str(length(duration))])
%%
mdl = fitlm(duration,amplitude);
plot(mdl),title(['R^2: ' num2str(mdl.Rsquared.Ordinary)])
