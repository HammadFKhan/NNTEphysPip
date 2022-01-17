function plotTF(TimeFreq,LFP)
oscillators = TimeFreq.tf.oscillators;
% wavelet = TimeFreq.tf.wavelet;
tf = TimeFreq.tf;
interpLevel = 4;
% imagesc(-1000:1000,oscillators.tb.f,test'),colormap(jet),axis xy, colorbar
figure,imagesc(-1000:1000,oscillators.tb.f,interp2(mean(oscillators.tb.C,3),interpLevel)),caxis([0 0.5]),title('\beta-\theta Phase Coupling'),colormap(jet),axis xy,set(gcf, 'Position',  [100, 100, 500, 500])
figure,imagesc(-1000:1000,oscillators.tg.f,interp2(mean(oscillators.tg.C,3),interpLevel)),caxis([0 0.5]),title('\theta-\gamma Phase Coupling'),colormap(jet),axis xy,set(gcf, 'Position',  [100, 100, 500, 500])
figure,imagesc(-1000:1000,oscillators.bg.f,interp2(mean(oscillators.bg.C,3),interpLevel)),caxis([0 0.5]),title('\beta-\gamma Phase Coupling'),colormap(jet),axis xy,set(gcf, 'Position',  [100, 100, 500, 500])
% figure()
% mtheta = TimeFreq.mtheta*(wavelet.theta_f(end)-wavelet.theta_f(1))+wavelet.theta_f(1);
% subplot(1,3,1),imagesc(-1000:1000,wavelet.theta_f,interp2(mean(abs(wavelet.theta_cfs),3),interpLevel)),title('\theta Spectrogram'),colorbar,colormap(jet),axis xy
% subplot(1,3,2),imagesc(-1000:1000,wavelet.beta_f,interp2(mean(abs(wavelet.beta_cfs),3),interpLevel)),title('\beta Spectrogram'),colorbar,colormap(jet),axis xy
% subplot(1,3,3),imagesc(-1000:1000,wavelet.gamma_f,interp2(mean(abs(wavelet.gamma_cfs),3),interpLevel)),title('\gamma Spectrogram'),colorbar,colormap(jet),axis xy
%% trial-by-trail mapping
% for trial = 1:size(wavelet.theta_cfs,3)
%     thetaTrial(trial,:) = mean(abs(wavelet.theta_cfs(:,:,trial)),2);
%     betaTrial(trial,:) = mean(abs(wavelet.beta_cfs(:,:,trial)),2);
%     gammaTrial(trial,:) = mean(abs(wavelet.gamma_cfs(:,:,trial)),2);
% end
% thetaPower = mean(thetaTrial)';
% betaPower = mean(betaTrial)';
% gammaPower = mean(gammaTrial)';
% figure,
% subplot(1,3,1),imagesc(thetaTrial),colormap(jet),colorbar, axis xy,caxis([0 25])
% subplot(1,3,2),imagesc(betaTrial),colormap(jet),colorbar, axis xy,caxis([0 25])
% subplot(1,3,3),imagesc((gammaTrial)),colormap(jet),colorbar, axis xy,caxis([0 25])


figure(),
% plotpad = size(tf.theta.theta,2);
% polarplot([zeros(1,plotpad); tf.theta.theta],[zeros(1,plotpad); ones(1,plotpad)],'k');
subplot(3,1,1),polarhistogram(horzcat(tf.theta.theta{:}));
title(['Theta ITPC: ' num2str(mean(tf.theta.itpc))]);
subplot(3,1,2),polarhistogram(horzcat(tf.beta.beta{:}),20);
title(['Beta ITPC: ' num2str(mean(tf.beta.itpc))]);
subplot(3,1,3),polarhistogram(horzcat(tf.gamma.gamma{:}),20);
title(['Gamma ITPC: ' num2str(mean(tf.gamma.itpc))]);


figure,
itpcAll = vertcat(tf.theta.itpc,tf.beta.itpc,tf.gamma.itpc);
itpcF = vertcat(tf.theta.f',tf.beta.f',tf.gamma.f');
plot(itpcF,itpcAll);