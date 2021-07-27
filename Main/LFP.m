%% plot LFP
function LFP(dataLow,channelStart,channelEnd)
fn = fieldnames(dataLow);
if nargin < 3, channelEnd = length(fn); end
if nargin < 2, channelStart = 1; end
spacing = 5000;
ycenter =[spacing:spacing:spacing*length(fn)]; 
colorgrad = zeros(length(fn),3);
colorgrad(:,1) = linspace(0.6350,0     ,length(fn));
colorgrad(:,2) = linspace(0.0780,0.4470,length(fn));
colorgrad(:,3) = linspace(0.1840,0.7410,length(fn));

for i = channelStart:channelEnd
    %plot LFP
    timelabel = 0:1/20000:60-1/20000;
    LFPData = dataLow.(string(fn(i)));
    plot(timelabel,LFPData'+ycenter(i),'color',colorgrad(i,:)), hold on
%     plot(mean(LFPData)+ycenter(i),'color','k','LineWidth',2)
    
%     %slightly under it, plot MUA
%     channelData = dataHigh.(['Channel' num2str(cc)]);
%     plot(trialTime,channelData'+ycenter(cc)-200,'color',colorgrad(cc,:))
%     plot(trialTime,mean(channelData)+ycenter(cc)-200,'color','k','LineWidth',2)    
end
% 
% rectangle('Position',[windowBefore./Fs.*1000,spacing,length(poleMoveIndices)/Fs*1000,spacing*length(channels)],...
%           'FaceColor',[0.25 0.25 0.25 0.2],...
%           'EdgeColor','none')
% alpha(0.2)
axis tight
box off 
set(gca,'YTick','','YTickLabel','')
title('LFP Freq Band')
xlabel('Time (s)')
ylabel('Channels')
disp('Done!')
