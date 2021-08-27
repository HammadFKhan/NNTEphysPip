function plotTF(Spikes,LFP)
tf = Spikes.tf;

figure,
subplot(4,1,1),stack_plot(LFP.LFP(:,1:2048)); title('Microprobe Data')
subplot(4,1,2),plot_matrix(tf.mtheta,tf.t,tf.theta.f); colormap(jet), title('Theta Spectrogram')
subplot(4,1,3),plot_matrix(tf.mbeta,tf.t,tf.beta.f); title('Beta Spectrogram')
subplot(4,1,4),plot_matrix(tf.mgamma,tf.t,tf.gamma.f); title('Gamma Spectrogram')

figure,
subplot(3,1,1),imagesc(0:tf.t,tf.theta.f,tf.theta.C'); clim([0,1]),ylabel('Frequency (Hz)'),axis xy,title('Theta Band Coherency'),colorbar,colormap(jet);
subplot(3,1,2),imagesc(0:tf.t,tf.beta.f,tf.beta.C'),ylabel('Frequency (Hz)')...
    ,axis xy,title('Beta Band Coherency'),colorbar, clim([0 1]);
subplot(3,1,3),imagesc(0:tf.t,tf.gamma.f,tf.gamma.C'),ylabel('Frequency (Hz)')...
    ,axis xy,title('Gamma Band Coherency'),colorbar, clim([0 1]);

figure,
% plotpad = size(tf.theta.theta,2);
% polarplot([zeros(1,plotpad); tf.theta.theta],[zeros(1,plotpad); ones(1,plotpad)],'k');
subplot(1,3,1),polarhistogram(tf.theta.theta,20);
title(['Theta ITPC: ' num2str(tf.theta.itpc)]);
subplot(1,3,2),polarhistogram(tf.beta.beta,20);
title(['Beta ITPC: ' num2str(tf.beta.itpc)]);
subplot(1,3,3),polarhistogram(tf.gamma.gamma,20);
title(['Gamma ITPC: ' num2str(tf.gamma.itpc)]);