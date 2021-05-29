%% plot MUA
function MUA(dataHigh,channelStart,channelEnd)
fn = fieldnames(dataHigh);
if nargin < 3, channelEnd = length(fn); end
if nargin < 2, channelStart = 1; end
spacing = 100;
ycenter =[spacing:spacing:spacing*length(fn)]; 
colorgrad = zeros(length(fn),3);
colorgrad(:,1) = linspace(0.6350,0     ,length(fn));
colorgrad(:,2) = linspace(0.0780,0.4470,length(fn));
colorgrad(:,3) = linspace(0.1840,0.7410,length(fn));

for i = channelStart:channelEnd
    channelData = dataHigh.(string(fn(i)));
    plot(channelData'+ycenter(i),'color',colorgrad(i,:)), hold on
    plot(mean(channelData)+ycenter(i),'color','k','LineWidth',2)
%         if(tt==numTrials)
%             plot(poleTime,ycenter(cc).*ones(1,length(poleTime)),'color',[0.25, 0.25, 0.25],'LineWidth',5)
%         end
end

% rectangle('Position',[windowBefore./Fs.*1000,spacing,length(poleMoveIndices)/Fs*1000,spacing*length(channels)],...
%           'FaceColor',[0.25 0.25 0.25 0.2],...
%           'EdgeColor','none')
% alpha(0.2)

axis tight
box off 
set(gca,'YTick','','YTickLabel','')
title('MUA Activity')
xlabel('Time (ms)')
ylabel('Channels')
disp('Done!')