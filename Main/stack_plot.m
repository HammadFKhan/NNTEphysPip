function stack_plot(DeltaFoverF)
Fs = 1024;

x = length(DeltaFoverF(1,:));
y = length(DeltaFoverF(:,1));
baseline = max(DeltaFoverF,[],'all');
time = (1:x)/Fs;
for i = 1:y
    gradient = i/y;
    plot(time,1*DeltaFoverF(i,:)+(.7*baseline),'LineWidth',1,'Color',[.5 .5 .5 .8]); hold on;
    baseline = baseline + max(DeltaFoverF,[],'all');
end
axis tight,box off
% set(gca,'XTick',[])
set(gca,'YTick',[])
disp('Done!')
end
