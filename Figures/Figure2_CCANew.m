%% CCA Conventional Analysis
mod1 = [];
mod2 = [];
mod3 = [];
for n = 1:6
    CC.hit.mod1(:,:,n) = combinedCCA(n).CCA.cueHit.CorrMap(:,:,1);
    CC.hit.mod2(:,:,n) = combinedCCA(n).CCA.cueHit.CorrMap(:,:,2);
    CC.hit.mod3(:,:,n) = combinedCCA(n).CCA.cueHit.CorrMap(:,:,3);
    
    CC.miss.mod1(:,:,n) = combinedCCA(n).CCA.cueMiss.CorrMap(:,:,1);
    CC.miss.mod2(:,:,n) = combinedCCA(n).CCA.cueMiss.CorrMap(:,:,2);
    CC.miss.mod3(:,:,n) = combinedCCA(n).CCA.cueMiss.CorrMap(:,:,3);
    
    CC.MIHit.mod1(:,:,n) = combinedCCA(n).CCA.MIHit.CorrMap(:,:,1);
    CC.MIHit.mod2(:,:,n) = combinedCCA(n).CCA.MIHit.CorrMap(:,:,2);
    CC.MIHit.mod3(:,:,n) = combinedCCA(n).CCA.MIHit.CorrMap(:,:,3);
    
    CC.MIFA.mod1(:,:,n) = combinedCCA(n).CCA.MIFA.CorrMap(:,:,1);
    CC.MIFA.mod2(:,:,n) = combinedCCA(n).CCA.MIFA.CorrMap(:,:,2);
    CC.MIFA.mod3(:,:,n) = combinedCCA(n).CCA.MIFA.CorrMap(:,:,3);
end

%%
sz = 6;
colors = [0 0.4470 0.7410;0.75 0.75 0.75;190/255 30/255 45/255];
time = linspace(-1.5,1.5,147);
f = figure;
f.Position = [681 559 960 420];
subplot(131),
dat = smoothdata(squeeze(mean(CC.hit.mod1(51,:,:),3)),'gaussian',10);
err = smoothdata(squeeze(std(CC.hit.mod1(50,:,:),[],3)),'gaussian',10);


plot(time,dat,'color',colors(1,:)),hold on
plot(time,dat+err/sqrt(sz),'color',colors(1,:))
plot(time,dat-err/sqrt(sz),'color',colors(1,:))


dat = smoothdata(squeeze(mean(CC.miss.mod1(51,:,:),3)),'gaussian',10);
err = smoothdata(squeeze(std(CC.miss.mod1(50,:,:),[],3)),'gaussian',10);
plot(time,dat,'color',colors(2,:)),hold on
plot(time,dat+err/sqrt(sz),'color',colors(2,:))
plot(time,dat-err/sqrt(sz),'color',colors(2,:))

box off,set(gca,'tickdir','out','fontsize',12),axis square,xlim([-0.5 1.5]),ylim([0.11 0.31])
ylabel('CC Coefficient')


subplot(132),
dat = smoothdata(squeeze(mean(CC.hit.mod2(51,:,:),3)),'gaussian',10);
err = smoothdata(squeeze(std(CC.hit.mod2(51,:,:),[],3)),'gaussian',10);


plot(time,dat,'color',colors(1,:)),hold on
plot(time,dat+err/sqrt(sz),'color',colors(1,:))
plot(time,dat-err/sqrt(sz),'color',colors(1,:))

dat = smoothdata(squeeze(mean(CC.miss.mod2(51,:,:),3)),'gaussian',10);
err = smoothdata(squeeze(std(CC.miss.mod2(51,:,:),[],3)),'gaussian',10);
plot(time,dat,'color',colors(2,:)),hold on
plot(time,dat+err/sqrt(sz),'color',colors(2,:))
plot(time,dat-err/sqrt(sz),'color',colors(2,:))

box off,set(gca,'tickdir','out','fontsize',12),axis square,xlim([-0.5 1.5]),ylim([0.11 0.31])
ylabel('CC Coefficient')


subplot(133),
dat = smoothdata(squeeze(mean(CC.hit.mod3(51,:,:),3)),'gaussian',10);
err = smoothdata(squeeze(std(CC.hit.mod3(51,:,:),[],3)),'gaussian',10);


plot(time,dat,'color',colors(1,:)),hold on
plot(time,dat+err/sqrt(sz),'color',colors(1,:))
plot(time,dat-err/sqrt(sz),'color',colors(1,:))

dat = smoothdata(squeeze(mean(CC.miss.mod3(51,:,:),3)),'gaussian',10);
err = smoothdata(squeeze(std(CC.miss.mod3(51,:,:),[],3)),'gaussian',10);
plot(time,dat,'color',colors(2,:)),hold on
plot(time,dat+err/sqrt(sz),'color',colors(2,:))
plot(time,dat-err/sqrt(sz),'color',colors(2,:))

box off,set(gca,'tickdir','out','fontsize',12),axis square,xlim([-0.5 1.5]),ylim([0.11 0.31])
ylabel('CC Coefficient')

%%% HIT AND FA

f = figure;
f.Position = [681 559 960 420];
subplot(131),
dat = smoothdata(squeeze(mean(CC.MIHit.mod1(51,:,:),3)),'gaussian',10);
err = smoothdata(squeeze(std(CC.MIHit.mod1(50,:,:),[],3)),'gaussian',10);


plot(time,dat,'color',colors(1,:)),hold on
plot(time,dat+err/sqrt(sz),'color',colors(1,:))
plot(time,dat-err/sqrt(sz),'color',colors(1,:))


dat = smoothdata(squeeze(mean(CC.MIFA.mod1(51,:,:),3)),'gaussian',10);
err = smoothdata(squeeze(std(CC.MIFA.mod1(50,:,:),[],3)),'gaussian',10);
plot(time,dat,'color',colors(3,:)),hold on
plot(time,dat+err/sqrt(sz),'color',colors(3,:))
plot(time,dat-err/sqrt(sz),'color',colors(3,:))

box off,set(gca,'tickdir','out','fontsize',12),axis square,xlim([-0.5 1.5]),ylim([0.155 0.31])
ylabel('CC Coefficient')


subplot(132),
dat = smoothdata(squeeze(mean(CC.MIHit.mod2(51,:,:),3)),'gaussian',10);
err = smoothdata(squeeze(std(CC.MIHit.mod2(51,:,:),[],3)),'gaussian',10);


plot(time,dat,'color',colors(1,:)),hold on
plot(time,dat+err/sqrt(sz),'color',colors(1,:))
plot(time,dat-err/sqrt(sz),'color',colors(1,:))

dat = smoothdata(squeeze(mean(CC.MIFA.mod2(51,:,:),3)),'gaussian',10);
err = smoothdata(squeeze(std(CC.MIFA.mod2(51,:,:),[],3)),'gaussian',10);
plot(time,dat,'color',colors(3,:)),hold on
plot(time,dat+err/sqrt(sz),'color',colors(3,:))
plot(time,dat-err/sqrt(sz),'color',colors(3,:))

box off,set(gca,'tickdir','out','fontsize',12),axis square,xlim([-0.5 1.5]),ylim([0.155 0.31])
ylabel('CC Coefficient')


subplot(133),
dat = smoothdata(squeeze(mean(CC.MIHit.mod3(51,:,:),3)),'gaussian',10);
err = smoothdata(squeeze(std(CC.MIHit.mod3(51,:,:),[],3)),'gaussian',10);


plot(time,dat,'color',colors(1,:)),hold on
plot(time,dat+err/sqrt(sz),'color',colors(1,:))
plot(time,dat-err/sqrt(sz),'color',colors(1,:))

dat = smoothdata(squeeze(mean(CC.MIFA.mod3(51,:,:),3)),'gaussian',10);
err = smoothdata(squeeze(std(CC.MIFA.mod3(51,:,:),[],3)),'gaussian',10);
plot(time,dat,'color',colors(3,:)),hold on
plot(time,dat+err/sqrt(sz),'color',colors(3,:))
plot(time,dat-err/sqrt(sz),'color',colors(3,:))

box off,set(gca,'tickdir','out','fontsize',12),axis square,xlim([-0.5 1.5]),ylim([0.155 0.31])
ylabel('CC Coefficient')
%%
argIn = combinedCCA(1).argIn;
mapDim = size(combinedCCA(1).CCA.cueHit.CorrMap, 2);
delays = (-argIn.MaxDelay:argIn.MaxDelay); % Convert to ms
t = (argIn.WindowLength/2:argIn.TimeStep:argIn.WindowLength/2+argIn.TimeStep*(mapDim-1)); % 

CANONICAL_PAIR_IDX = 1;
baselinePeriodStart = 500;
evokedPeriodStart = 1500;
evokedPeriodStop = 2500;
baselineIdx = find(t>=baselinePeriodStart,1);
startIdx = find(t>=evokedPeriodStart,1);
stopIdx = find(t>=evokedPeriodStop,1);

delayCCA.cueHit = zeros(size(delays,2),size(combinedCCA,2));

for i=1:size(delays,2)
    delayCCA.cueHit(i,:) = cell2mat(arrayfun(@(s) mean(s.CCA.cueHit.CorrMap(i,startIdx:stopIdx,CANONICAL_PAIR_IDX))-mean(s.CCA.cueHit.CorrMap(baselineIdx,1:startIdx-1,CANONICAL_PAIR_IDX)), combinedCCA,'UniformOutput',false));
end

figure();
dat = smoothdata(squeeze(mean(delayCCA.cueHit,2)),'gaussian',9);
err = smoothdata((squeeze(std(delayCCA.cueHit,1,2)))/sqrt(size(delayCCA.cueHit,2)),'gaussian',9);
h1=plot(delays,dat,'Color', [0 0.1 0.8],'LineWidth',2); hold on;
plot(delays,dat+err ,'Color', [0 0.1 0.8 0.4],'LineWidth',2);
plot(delays,dat-err ,'Color', [0 0.1 0.8 0.4],'LineWidth',2);
xlabel('Delay (ms)'); ylabel('Evoked CCA');
title('Evoked CCA wrt to delays')
box off;set(gca,'TickDir','out','fontsize',14');
drawnow;


