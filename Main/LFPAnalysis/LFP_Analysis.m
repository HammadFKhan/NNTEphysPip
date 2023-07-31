function filtData = LFP_Analysis(filtData,Ripples)
Fs = 20000;
LFP = filtData.lowpassData(filtData.bestLFPchan,:);
LFPm = sgolayfilt(LFP,9,2999);

figure,plot(LFP(40000:80000)); hold on;
plot(LFPm(40000:80000),'LineWidth',4);ylim([-500 500])
