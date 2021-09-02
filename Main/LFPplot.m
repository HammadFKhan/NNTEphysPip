function LFPplot(LFP)
Fs = LFP.downSampleFreq;
win = 1;
figure(3),
 %Plot the first 5 seconds of data
subplot(4,1,1),plot(LFP.times(:,1:(win*Fs)),LFP.bestLFP(:,1:(win*Fs)),'k'); title('Raw Data'),box off
subplot(4,1,2),plot(LFP.times(:,1:(win*Fs)),LFP.theta_band(:,1:(win*Fs)),'r');title('\theta Band'),box off
subplot(4,1,3),plot(LFP.times(:,1:(win*Fs)),LFP.beta_band(:,1:(win*Fs)),'m');title('\beta Band'),box off
subplot(4,1,4),plot(LFP.times(:,1:(win*Fs)),LFP.gamma_band(:,1:(win*Fs)),'g');title('\gamma Band'),box off
xlabel('Time (s)'), ylabel('Voltage (\muV)')

figure(4)
subplot(3,1,1),plot(LFP.times(:,1:(win*Fs)),LFP.theta_temppow(:,1:(win*Fs)),'r');title('\theta Band'),box off
subplot(3,1,2),plot(LFP.times(:,1:(win*Fs)),LFP.beta_temppow(:,1:(win*Fs)),'m');title('\beta Band'),box off
subplot(3,1,3),plot(LFP.times(:,1:(win*Fs)),LFP.beta_temppow(:,1:(win*Fs)),'g');title('\gamma Band'),box off
xlabel('Time (s)'), ylabel('Power (dB)')
