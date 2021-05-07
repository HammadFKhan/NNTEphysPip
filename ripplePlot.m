function ripplePlot(rippleData,data,stats,maps)
% Optionally, plot results
figure('Name','Detected Signal')
plot(rippleData.timestamps,rippleData.signal);hold on;
for j=1:size(rippleData.ripples,1)
    plot([rippleData.ripples(j,1) rippleData.ripples(j,1)],ylim,'g-');
    plot([rippleData.ripples(j,2) rippleData.ripples(j,2)],ylim,'k-');
    plot([rippleData.ripples(j,3) rippleData.ripples(j,3)],ylim,'r-');
end
plot(xlim,[rippleData.lowThresholdFactor rippleData.lowThresholdFactor],'k','linestyle','--');
plot(xlim,[rippleData.highThresholdFactor rippleData.highThresholdFactor],'k-');
axis tight
box off
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
    print(gcf,'-painters','-depsc', 'Figures/RippleDetect.eps', '-r250');

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
plot(((1:rippleData.nBins)'-ceil(rippleData.nBins/2))/rippleData.nBins*diff(rippleData.durations),maps(electrode).ripples','b');

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