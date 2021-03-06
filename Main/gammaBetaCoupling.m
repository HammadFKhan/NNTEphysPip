function betaGammaCoupling = gammaBetaCoupling(mGamma,mBeta)
%% Gamma amplitude phase coupling by depth and trial
if size(mGamma) & size(mBeta) > 2 % Beta coupling by depth
    for ii = 1:size(mBeta,2)
%         figure,hold on
        for j = 1:size(mBeta,3)
        mGammaAmplitude = mGamma(:,ii,j).^2;
        mBetaPhase = angle(hilbert(mBeta(:,ii,j)));
        %discretize the angle into 18 equal bins
        dBetaPhase = discretize(mBetaPhase*180/pi,-180:20:180);
        % Sort gamma into beta bins
        for i = 1:max(dBetaPhase)
            betaGammaCoupling(i,ii,j) = median(mGammaAmplitude(dBetaPhase==i,1));
        end
%         bar(-179:20:180,betaGammaCoupling)
%         xlim([-180 180])
        end
    end
else
    for ii = 1:size(mBeta,2)
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
end
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