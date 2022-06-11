%IEI stats
%Percentage state
data = [ -diff(l23rest) -diff(l23init) -diff(l23run);-diff(l5rest) -diff(l5init) -diff(l5run)]; %we use -diff because it is burst-single
x = [23 5];
figure,bar(x, data),title('W12'),ylim([-1 1])
legend('Init','Run','Rest')
%%
data = L23RestEventISI(L23RestEventISI>0);
figure,boxplot(data,'PlotStyle','compact'),hold on,title('L23Rest'),box off
plot(1.2*ones(length(data),1),data,'.'),ylim([0 0.9])

data = L23InitiateEventISI(L23InitiateEventISI>0);
figure,boxplot(data,'PlotStyle','compact'),hold on,title('L23Initiate'),box off
plot(1.2*ones(length(data),1),data,'.'),ylim([0 0.9])

data = L23RunEventISI(L23RunEventISI>0);
figure,boxplot(data,'PlotStyle','compact'),hold on,title('L23Run'),box off
plot(1.2*ones(length(data),1),data,'.'),ylim([0 0.9])

data = L5RestEventISI(L5RestEventISI>0);
figure,boxplot(data,'PlotStyle','compact'),hold on,title('L5Rest'),box off
plot(1.2*ones(length(data),1),data,'.'),ylim([0 0.9])

data = L5InitiateEventISI(L5InitiateEventISI>0);
figure,boxplot(data,'PlotStyle','compact'),hold on,title('L5Initiate'),box off
plot(1.2*ones(length(data),1),data,'.'),ylim([0 0.9])

data = L5RunEventISI(L5RunEventISI>0);
figure,boxplot(data,'PlotStyle','compact'),hold on,title('L5Run'),box off
plot(1.2*ones(length(data),1),data,'.'),ylim([0 0.9])