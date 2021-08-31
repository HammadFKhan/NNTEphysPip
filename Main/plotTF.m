function plotTF(Spikes,LFP)
tf = Spikes.tf;

% figure,plot_matrix(tf.mbroad,tf.t,tf.broad.f); colormap(jet), title('Theta Spectrogram')
% figure,imagesc(0:tf.t,tf.broad.f,tf.broad.C(:,:,100)'),ylabel('Frequency (Hz)')...
%     ,axis xy,title('Beta Band Coherency'),colorbar;
% figure,imagesc(0:tf.t,tf.broad.f,mean(tf.broad.phi,3)'),ylabel('Frequency (Hz)')...
%     ,axis xy,title('Beta Band Phase'),colorbar;
figure,
subplot(4,1,1),stack_plot(LFP.LFP(1:12,1:2048)); title('Microprobe Data')
subplot(4,1,2),plot_matrix(tf.mtheta,tf.t,tf.theta.f); colormap(jet), title('Theta Spectrogram')
subplot(4,1,3),plot_matrix(tf.mbeta,tf.t,tf.beta.f); title('Beta Spectrogram')
subplot(4,1,4),plot_matrix(tf.mgamma,tf.t,tf.gamma.f); title('Gamma Spectrogram')

figure,
interpLevel = 3;
subplot(3,1,1),imagesc(-1000:1000,tf.theta.f,interp2(tf.theta.C(:,:,1),interpLevel)');ylabel('Frequency (Hz)'),...
    axis xy,title('Theta Band Coherency'),colorbar,colormap(jet);
subplot(3,1,2),imagesc(0:tf.t,tf.beta.f,interp2(mean(tf.beta.C,3),interpLevel)'),ylel('Frequency (Hz)')...
    ,axis xy,title('Beta Band Coherency'),colorbar;
subplot(3,1,3),imagesc(0:tf.t,tf.gamma.f,interp2(mean(tf.gamma.C,3),interpLevel)'),ylabel('Frequency (Hz)')...
    ,axis xy,title('Gamma Band Coherency'),colorbar;

figure,
subplot(3,1,1),imagesc(0:tf.t,tf.theta.f,mean(tf.theta.phi,3)');ylabel('Frequency (Hz)'),axis xy,title('Theta Band Phase'),colorbar,colormap(jet);
subplot(3,1,2),imagesc(0:tf.t,tf.beta.f,mean(tf.beta.phi,3)'),ylabel('Frequency (Hz)')...
    ,axis xy,title('Beta Band Phase'),colorbar;
subplot(3,1,3),imagesc(0:tf.t,tf.gamma.f,mean(tf.gamma.phi,3)'),ylabel('Frequency (Hz)')...
    ,axis xy,title('Gamma Band Phase'),colorbar;
figure,
imagesc(-1000:1000,tf.wavelet.f,mean(abs(tf.wavelet.cfs),3)),title('Wavelet Spectrogram'),colormap(jet),axis xy
yline(2,'w--','Linewidth',2)
yline(10,'w--','Linewidth',2)
yline(30,'w--','Linewidth',2)
yline(80,'w--','Linewidth',2)
figure,
subplot(3,1,1),imagesc(-1000:1000,tf.wavelet.theta_f,mean(abs(tf.wavelet.theta_cfs),3)),title('Wavelet Spectrogram'),colorbar,colormap(jet),axis xy
subplot(3,1,2),imagesc(-1000:1000,tf.wavelet.beta_f,mean(abs(tf.wavelet.beta_cfs),3)),title('Wavelet Spectrogram'),colorbar,colormap(jet),axis xy
subplot(3,1,3),imagesc(-1000:1000,tf.wavelet.gamma_f,mean(abs(tf.wavelet.gamma_cfs),3)),title('Wavelet Spectrogram'),colorbar,colormap(jet),axis xy
% figure,
% % plotpad = size(tf.theta.theta,2);
% % polarplot([zeros(1,plotpad); tf.theta.theta],[zeros(1,plotpad); ones(1,plotpad)],'k');
% subplot(1,3,1),polarhistogram(tf.theta.theta,20);
% title(['Theta ITPC: ' num2str(tf.theta.itpc)]);
% subplot(1,3,2),polarhistogram(tf.beta.beta,20);
% title(['Beta ITPC: ' num2str(tf.beta.itpc)]);
% subplot(1,3,3),polarhistogram(tf.gamma.gamma,20);
% title(['Gamma ITPC: ' num2str(tf.gamma.itpc)]);