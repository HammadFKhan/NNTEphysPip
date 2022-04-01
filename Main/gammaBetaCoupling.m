function betaGammaCoupling = gammaBetaCoupling(mGamma,mBeta)
%% Gamma amplitude phase coupling
figure,hold on
for ii = 1:11
mGammaAmplitude = mGamma(:,ii).^2;
mBetaPhase = angle(hilbert(mBeta(:,ii)));
%discretize the angle into 18 equal bins
dBetaPhase = discretize(mBetaPhase*180/pi,-180:20:180);
% Sort gamma into beta bins
for i = 1:max(dBetaPhase)
    betaGammaCoupling(i) = mean(mGammaAmplitude(dBetaPhase==i,1));
end
bar(-179:20:180,betaGammaCoupling)
end
xlim([-180 180])
%% Random Plots

% data = mGamma(:,1);
% figure,
% subplot(311)
% plot(data)
% xlabel('Time (ms)'), ylabel('Amplitude (\muV)')
% 
% % plot power (squared magnitude from origin to dot-product location in
% % complex space)
% subplot(312)
% plot(data.^2)
% xlabel('Time (ms)'), ylabel('Power \muV^2')
% 
% 
% % plot phase (angle of vector to dot-product, relative to positive real
% % axis)
% subplot(313)
% plot(angle(hilbert(data)))
% xlabel('Time (ms)'), ylabel('Phase (rad.)')
% 
% data = mBeta(:,1);
% figure,
% subplot(311)
% plot(data)
% xlabel('Time (ms)'), ylabel('Amplitude (\muV)')
% 
% % plot power (squared magnitude from origin to dot-product location in
% % complex space)
% subplot(312)
% plot(data.^2)
% xlabel('Time (ms)'), ylabel('Power \muV^2')
% 
% 
% % plot phase (angle of vector to dot-product, relative to positive real
% % axis)
% subplot(313)
% plot(angle(hilbert(data)))
% xlabel('Time (ms)'), ylabel('Phase (rad.)')