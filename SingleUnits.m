%% plot single unit waveforms
function SingleUnits(allSpikes,Fs,spacing,channelStart,channelEnd)
fn = fieldnames(allSpikes);
if nargin < 5, channelEnd = length(fn); end
if nargin < 4, channelStart = 1; end
if nargin < 3, spacing = 100; end
if nargin < 2, Fs = 20000; end
space = 0;
for i = channelStart:channelEnd
    theseSpikes = allSpikes.(string(fn(i)));
    ycenter =[spacing:spacing:spacing*length(theseSpikes)]; 
    t = (1:length(theseSpikes(1,:)))./Fs.*1000;
    plot(t,theseSpikes+space,'Color',[0.5 0.5 0.5 0.4]), hold on
    plot(t,mean(theseSpikes)+space,'k','LineWidth',2)    
    space = space+spacing;
    axis tight
    box off
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    
end
disp('Done!')
