function stack_plot(DeltaFoverF,space,scale)
if nargin<3,scale=1; end
if nargin<2,space=1; end

Fs = 1024;

x = length(DeltaFoverF(1,:));
y = length(DeltaFoverF(:,1));
baseline = max(DeltaFoverF,[],'all');
time = (1:x)/Fs;
for i = 1:y
    gradient = i/y;
    plot(time,scale*DeltaFoverF(i,:)+(space*baseline),'LineWidth',1,'Color',[.5 .5 .5 .8]); hold on;
    baseline = baseline + max(DeltaFoverF,[],'all');
end
axis tight,box off
% set(gca,'XTick',[])
set(gca,'YTick',[])
disp('Done!')
end
