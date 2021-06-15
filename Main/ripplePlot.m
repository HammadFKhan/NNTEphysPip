function ripplePlot(Ripples)
% Parsing
maps = Ripples.maps;
stats = Ripples.stats;
data = Ripples.data;
LFP = Ripples.rippleOnset.LFP;
SWR = Ripples.rippleOnset.SWR;
cfs = Ripples.rippleOnset.cfs;
waveletFreq = Ripples.rippleOnset.f;
% Optionally, plot results
% figure('Name','Detected Signal')
% subplot(4,1,[1 3]),plot(Ripples.timestamps,Ripples.signal);hold on;
% for j=1:size(Ripples.ripples,1)
%     plot([Ripples.ripples(j,1) Ripples.ripples(j,1)],ylim,'g-');
%     plot([Ripples.ripples(j,2) Ripples.ripples(j,2)],ylim,'k-');
%     plot([Ripples.ripples(j,3) Ripples.ripples(j,3)],ylim,'r-');
% end
% plot(xlim,[Ripples.lowThresholdFactor Ripples.lowThresholdFactor],'k','linestyle','--');
% plot(xlim,[Ripples.highThresholdFactor Ripples.highThresholdFactor],'k-');
% subplot(4,1,4),plot(VR_data.AvgVel{1,1});
% axis tight
% box off
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
%     print(gcf,'-painters','-depsc', 'Figures/RippleDetect.eps', '-r250');

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
    print(gcf,'-painters','-dpng', 'Figures/Ripple.png', '-r250');

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
plot(((1:Ripples.nBins)'-ceil(Ripples.nBins/2))/Ripples.nBins*diff(Ripples.durations),maps(electrode).ripples','b');

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

%% Ripple Onset

figure('Name','SWR Onset Frequency')
for ii = 1:size(SWR,2)
    plotSize = ceil(sqrt(size(SWR,2)));
    subplot(plotSize,plotSize,ii),imagesc(-150:150,waveletFreq,abs(cfs(:,:,ii)))
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    axis xy
    colormap(jet)
    % ylim([0 250])
%     title('CWT of Ripple Data')
end
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
%     print(gcf,'-painters','-depsc', 'Figures/SWRonset.eps', '-r250');

% figure('Name','CSD');
% channelmap = ones(size(LFP,2),1);
% avgLFP = mean(LFP,3);
% try
%     Vq = interp2(avgLFP,5);
% catch ME
%     disp('Sample set is too small for interpolation, plotting raw...');
%     Vq = avgLFP;
% end
% imagesc(-150:150,channelmap,Vq);colormap(jet); colorbar;box off;set(gca,'YTick',[]);
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
%     print(gcf,'-painters','-depsc', 'Figures/CSD.eps', '-r250');
end