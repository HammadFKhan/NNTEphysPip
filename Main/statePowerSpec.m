%% Calculate broadband beta power
function statePowerSpec(TimeFreq,LFP)
beta = TimeFreq.beta;
theta = TimeFreq.theta;
gamma = TimeFreq.gamma;
betaM = squeeze(mean(beta,2));
thetaM = squeeze(mean(theta,2));
gammaM = squeeze(mean(gamma,2));
betaM = betaM';thetaM = thetaM';gammaM = gammaM';
refSignal = mean(LFP.LFP,1);
for i = 1:size(betaM,1)
    [pxxRef(:,i),f] = pwelch(refSignal((i-1)*size(betaM,2)+1:i*size(betaM,2)),500,300,500,1024);
    [pxx(:,i),~] = pwelch(betaM(i,:),500,300,500,1024);
    [pxx(:,i),~] = pwelch(betaM(i,:),500,300,500,1024);
    [pxx1(:,i),~] = pwelch(thetaM(i,:),500,300,500,1024);
    [pxx2(:,i),~] = pwelch(gammaM(i,:),500,300,500,1024);
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
figure,semilogx(f,10*log10(tpxx1.*tpxx2.*tpxx.*tpxxRef)),xlim([0 100])