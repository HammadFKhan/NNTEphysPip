function plotTF(TimeFreq,LFP)
oscillators = TimeFreq.tf.oscillators;
wavelet = TimeFreq.tf.wavelet;
interpLevel = 2;
figure(5)
% imagesc(-1000:1000,oscillators.tb.f,test'),colormap(jet),axis xy, colorbar
imagesc(-1000:1000,oscillators.tg.f,interp2(mean(oscillators.tb.phi,3)',interpLevel)),colormap(jet),axis xy, colorbar
imagesc(-1000:1000,oscillators.tg.f,interp2(mean(oscillators.tg.phi,3)',interpLevel)),colormap(jet),axis xy, colorbar
imagesc(-1000:1000,oscillators.tb.f,interp2(mean(oscillators.bg.phi,3)',interpLevel)),colormap(jet),axis xy, colorbar

figure(6)
mtheta = TimeFreq.mtheta*(wavelet.theta_f(end)-wavelet.theta_f(1))+wavelet.theta_f(1);
subplot(2,3,1),imagesc(-1000:1000,wavelet.theta_f,interp2(mean(abs(wavelet.theta_cfs),3),interpLevel)),title('\theta Spectrogram'),colormap(jet),axis xy
subplot(2,3,4),plot(TimeFreq.mtheta),axis off;
subplot(2,3,2),imagesc(-1000:1000,wavelet.beta_f,interp2(mean(abs(wavelet.beta_cfs),3),interpLevel)),title('\beta Spectrogram'),colormap(jet),axis xy
subplot(2,3,5),plot(TimeFreq.mbeta),axis off;
subplot(2,3,3),imagesc(-1000:1000,wavelet.gamma_f,interp2(mean(abs(wavelet.gamma_cfs),3),interpLevel)),title('\gamma Spectrogram'),colormap(jet),axis xy
subplot(2,3,6),plot(TimeFreq.mgamma),axis off;

% figure,
% % plotpad = size(tf.theta.theta,2);
% % polarplot([zeros(1,plotpad); tf.theta.theta],[zeros(1,plotpad); ones(1,plotpad)],'k');
% subplot(1,3,1),polarhistogram(tf.theta.theta,20);
% title(['Theta ITPC: ' num2str(tf.theta.itpc)]);
% subplot(1,3,2),polarhistogram(tf.beta.beta,20);
% title(['Beta ITPC: ' num2str(tf.beta.itpc)]);
% subplot(1,3,3),polarhistogram(tf.gamma.gamma,20);
% title(['Gamma ITPC: ' num2str(tf.gamma.itpc)]);

test = sort(oscillators.tb.phi,1);
test = test.*(180/pi);
 C = permute(test,[1 3 2]);
 C = reshape(C,[],size(test,2),1);
 edgesPhase = -180:1:180; %Bin width on the track
 phaseBin = discretize(C,edgesPhase);