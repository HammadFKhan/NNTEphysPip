%IEI stats
%Percentage state
data = [ -diff(l23rest) -diff(l23init) -diff(l23run);-diff(l5rest) -diff(l5init) -diff(l5run)]; %we use -diff because it is burst-single
x = [23 5];
figure,bar(x, data),title('W12'),ylim([-1 1])
legend('Init','Run','Rest')