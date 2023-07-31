%% Calculate broadband beta power
function pxSpecs = statePowerSpec(TimeFreq,LFP)
beta = TimeFreq.tfRun.betaLFP;
theta = TimeFreq.tfRun.thetaLFP;
gamma = TimeFreq.tfRun.gammaLFP;
full = TimeFreq.tfRun.fullLFP;
betaM = squeeze(mean(beta,2));
thetaM = squeeze(mean(theta,2));
gammaM = squeeze(mean(gamma,2));
fullM = squeeze(mean(full,2));
betaM = betaM';thetaM = thetaM';gammaM = gammaM';fullM = fullM';
refSignal = mean(LFP.LFP,1);
for i = 1:size(betaM,1)
    %     [pxxRef(:,i),f] = pwelch(refSignal((i-1)*size(betaM,2)+1:i*size(betaM,2)),500,300,1000,1024);
    [pxxRef(:,i),f] = pwelch(fullM(i,:),500,30,1000,1024);
    [pxx(:,i),~] = pwelch(betaM(i,:),500,30,1000,1024);
    [pxx1(:,i),~] = pwelch(thetaM(i,:),500,30,1000,1024);
    [pxx2(:,i),~] = pwelch(gammaM(i,:),500,30,1000,1024);
    tpxx(:,i) = pxx(:,i);%.*(f.^2);
    tpxx1(:,i) = pxx1(:,i);%.*(f1.^2);
    tpxx2(:,i) = pxx2(:,i);%.*(f2.^2);
    tpxxRef(:,i) = pxxRef(:,i).*f.^2;
end

% pxx(1,:)= 0;
figure,semilogx(f,mean(10*log10(tpxx),2)),hold on
semilogx(f,mean(10*log10(tpxx1),2))
semilogx(f,mean(10*log10(tpxx2),2))
semilogx(f,mean(10*log10(tpxxRef),2))
broadband = 10*log10(tpxx1.*tpxx2.*tpxx.*tpxxRef);
mbroadband = mean(broadband,2);
figure, hold on,semilogx(f,broadband),semilogx(f,mbroadband,'k','LineWidth',3),xlim([0 100]),set(gca, 'XScale', 'log')
figure,lineError(f,broadband'),set(gca, 'XScale', 'log'),xlim([0 100])



% broadband wavelet power spectrum
for i = 1:size(fullM,1)
% [waveletT,fT] = cwt(mean(thetaM,1),1024,'FrequencyLimit',[4 10]);
% [waveletB,fB] = cwt(mean(betaM,1),1024,'FrequencyLimit',[10 30]);
% [waveletG,fG] = cwt(mean(gammaM,1),1024,'FrequencyLimit',[30 80]);
[wavelet(:,:,i),fwav] = cwt(fullM(i,:),1024,'FrequencyLimit',[10 32]);
end

z = abs(mean(wavelet,3));
z = (z-mean(z,2))./(max(z,[],2)-min(z,[],2));
figure,imagesc(-1000:1000,fwav,abs(wavelet(:,:,7))),colormap(jet),axis xy,colorbar, hold on
figure,imagesc(-1000:1000,TimeFreq.tfRun.beta.f,TimeFreq.tfRun.depth.L23.beta.C(:,:,3)'),colormap(jet),axis xy,colorbar,hold on
L23Coherence = squeeze(mean(TimeFreq.tfRun.depth.L23.beta.C(15:end,:,:),1));
L23Coherence = squeeze(mean(L23Coherence,1));

L23Phi = squeeze(mean(TimeFreq.tfRun.depth.L23.beta.phi(15:end,:,1:10),1));
L23Phi = squeeze(mean(L23Phi,1));

L5Coherence = squeeze(mean(TimeFreq.tfRun.depth.L5.beta.C(15:end,:,:),1));
L5Coherence = squeeze(mean(L5Coherence,1));

L5Phi = squeeze(mean(TimeFreq.tfRun.depth.L5.beta.phi(15:end,:,1:10),1));
L5Phi = squeeze(mean(L5Phi,1));

figure,boxplot(L23Coherence,'Plotstyle','compact'),ylim([0 .7]),box off
figure,boxplot(L5Coherence,'Plotstyle','compact'),ylim([0 .7]),box off
figure,boxplot((180/pi)*L23Phi,'Plotstyle','compact'),ylim([-180 180]),box off
figure,boxplot((180/pi)*L5Phi,'Plotstyle','compact'),ylim([-180 180]),box off

%% Output 
pxSpecs.broadband = broadband;
pxSpecs.f = f;
pxSpecs.wavelet = wavelet;
pxSpecs.fwav = fwav;
pxSpecs.waveletNorm = z;


