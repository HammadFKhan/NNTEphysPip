function LFPplot(LFP)
Fs = LFP.downSampleFreq;
start = 200*Fs;
win = 2*Fs;
stop = start+win-1;
figure(3),set(gcf, 'Position',  [100, 100, 400, 500]),
 %Plot the first 5 seconds of data
subplot(4,1,1),plot(LFP.times(:,1:win),LFP.LFP(1,start:stop),'k'); title('Raw Data'),box off,ylim([-300 300])
subplot(4,1,2),plot(LFP.times(:,1:win),LFP.theta_band(:,start:stop),'r');title('\theta Band'),box off,ylim([-100 100])
subplot(4,1,3),plot(LFP.times(:,1:win),LFP.beta_band(:,start:stop),'m');title('\beta Band'),box off,ylim([-100 100])
subplot(4,1,4),plot(LFP.times(:,1:win),LFP.gamma_band(:,start:stop),'g');title('\gamma Band'),box off,ylim([-100 100])
xlabel('Time (s)'), ylabel('Voltage (\muV)')

figure(4)
subplot(3,1,1),plot(LFP.times(:,1:win),LFP.theta_temppow(:,start:stop),'r');title('\theta Band'),box off
subplot(3,1,2),plot(LFP.times(:,1:win),LFP.beta_temppow(:,start:(stop)),'m');title('\beta Band'),box off
subplot(3,1,3),plot(LFP.times(:,1:win),LFP.beta_temppow(:,start:stop),'g');title('\gamma Band'),box off
xlabel('Time (s)'), ylabel('Power (dB)')
