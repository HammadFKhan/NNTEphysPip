function plotTF(TimeFreq,LFP)
oscillators = TimeFreq.tf.oscillators;
wavelet = TimeFreq.tf.wavelet;
tf = TimeFreq.tf;
interpLevel = 2;
figure(),set(gcf, 'Position',  [100, 100, 500, 100]),
% imagesc(-1000:1000,oscillators.tb.f,test'),colormap(jet),axis xy, colorbar
subplot(1,3,1),imagesc(-1000:1000,oscillators.tg.f,interp2(mean(oscillators.tb.phi,3)',interpLevel)),colorbar,title('\theta Phase'),colormap(jet),axis xy,caxis([-0.35 0.35]) 
subplot(1,3,2),imagesc(-1000:1000,oscillators.tg.f,interp2(mean(oscillators.tg.phi,3)',interpLevel)),colorbar,title('\beta Phase'),colormap(jet),axis xy, caxis([-0.35 0.35])
subplot(1,3,3),imagesc(-1000:1000,oscillators.tb.f,interp2(mean(oscillators.bg.phi,3)',interpLevel)),colorbar,title('\gamma Phase'),colormap(jet),axis xy,caxis([-0.35 0.35])
ylim([10 30])
figure()
mtheta = TimeFreq.mtheta*(wavelet.theta_f(end)-wavelet.theta_f(1))+wavelet.theta_f(1);
subplot(1,3,1),imagesc(-1000:1000,wavelet.theta_f,interp2(mean(abs(wavelet.theta_cfs),3),interpLevel)),title('\theta Spectrogram'),colorbar,colormap(jet),axis xy
subplot(1,3,2),imagesc(-1000:1000,wavelet.beta_f,interp2(mean(abs(wavelet.beta_cfs),3),interpLevel)),title('\beta Spectrogram'),colorbar,colormap(jet),axis xy
subplot(1,3,3),imagesc(-1000:1000,wavelet.gamma_f,interp2(mean(abs(wavelet.gamma_cfs),3),interpLevel)),title('\gamma Spectrogram'),colorbar,colormap(jet),axis xy
%% trial-by-trail mapping
for trial = 1:size(wavelet.theta_cfs,3)
    thetaTrial(trial,:) = mean(abs(wavelet.theta_cfs(:,:,trial)),2);
    betaTrial(trial,:) = mean(abs(wavelet.beta_cfs(:,:,trial)),2);
    gammaTrial(trial,:) = mean(abs(wavelet.gamma_cfs(:,:,trial)),2);
end
thetaPower = mean(thetaTrial)';
betaPower = mean(betaTrial)';
gammaPower = mean(gammaTrial)';
figure,
subplot(1,3,1),imagesc(thetaTrial),colormap(jet),colorbar, axis xy,caxis([0 25])
subplot(1,3,2),imagesc(betaTrial),colormap(jet),colorbar, axis xy,caxis([0 25])
subplot(1,3,3),imagesc((gammaTrial)),colormap(jet),colorbar, axis xy,caxis([0 25])


% for i = 1:size(tf.theta.itpc,1)
% figure(10+i),
% % plotpad = size(tf.theta.theta,2);
% % polarplot([zeros(1,plotpad); tf.theta.theta],[zeros(1,plotpad); ones(1,plotpad)],'k');
% subplot(3,1,1),polarhistogram(tf.theta.theta{i},20);
% title(['Theta ITPC: ' num2str(tf.theta.itpc(i))]);
% subplot(3,1,2),polarhistogram(tf.beta.beta{i},20);
% title(['Beta ITPC: ' num2str(tf.beta.itpc(i))]);
% subplot(3,1,3),polarhistogram(tf.gamma.gamma{i},20);
% title(['Gamma ITPC: ' num2str(tf.gamma.itpc(i))]);
% end
% subplot(3,1,1),polarhistogram(tf.theta.theta,20);
% title(['Theta ITPC: ' num2str(tf.theta.itpc)]);
% subplot(3,1,2),polarhistogram(tf.beta.beta,20);
% title(['Beta ITPC: ' num2str(tf.beta.itpc)]);
% subplot(3,1,3),polarhistogram(tf.gamma.gamma,20);
% title(['Gamma ITPC: ' num2str(tf.gamma.itpc)]);
figure,
% subplot(1,3,1),plot(tf.theta.f',tf.theta.itpc)
% subplot(1,3,2),plot(tf.beta.f',tf.beta.itpc)
% subplot(1,3,3),plot(tf.gamma.f',tf.gamma.itpc)
figure,
itpcAll = vertcat(tf.theta.itpc,tf.beta.itpc,tf.gamma.itpc);
itpcF = vertcat(tf.theta.f',tf.beta.f',tf.gamma.f');
plot(itpcF,itpcAll);
% test = sort(oscillators.tb.phi,1);
% test = test.*(180/pi);
%  C = permute(test,[1 3 2]);
%  C = reshape(C,[],size(test,2),1);
%  edgesPhase = -180:1:180; %Bin width on the track
%  phaseBin = discretize(C,edgesPhase);