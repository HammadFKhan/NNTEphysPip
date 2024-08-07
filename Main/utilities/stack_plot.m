function stack_plot(DeltaFoverF,space,scale,Fs,normFlag)
if nargin < 5 || strcmp(normFlag,'')
    normFlag = 0;
end
if isempty(scale),scale=1; end
if isempty(space),space=1; end
if isempty(Fs),Fs = 1024; end

x = length(DeltaFoverF(1,:));
y = length(DeltaFoverF(:,1));
baseline = max(DeltaFoverF,[],'all');
time = (1:x)/Fs;
% [grad,~]=colorGradient([247 224 10]/255,[10 247 148]/255,y);
[grad,~]=colorGradient([7 49 97]/255,[110 192 235]/255,y);
DeltaFoverF = flip(DeltaFoverF,1); %Flip because of the way we plot it
for i = 1:y
    gradient = i/y;
    if normFlag
        DeltaFoverF(i,:) = (DeltaFoverF(i,:)-min(DeltaFoverF(i,:)))/(max(DeltaFoverF(i,:))-min(DeltaFoverF(i,:)));
    end
    plot(time,scale*DeltaFoverF(i,:)+(space*baseline),'LineWidth',1,'Color',grad(i,:)); hold on;
    baseline = baseline + max(DeltaFoverF,[],'all');
end
axis tight,box off
% set(gca,'XTick',[])
% set(gca,'YTick',[])
set(gca,'TickDir','out')
disp('Done!')
end
