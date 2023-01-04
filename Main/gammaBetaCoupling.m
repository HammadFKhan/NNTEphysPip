function betaGammaCoupling = gammaBetaCoupling(LFP,Tf,betaGroup)
mGamma = Tf.gammaLFP;
mBeta = Tf.betaLFP;
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
%% beta event/gamma amplitude modulation
% gamma amplitude during beta events
betaEvent = betaGroup(64).electrode.betaBurst.detectedBeta; %use global LFP values as beta bursts is recorded with this reference
beta = LFP.beta_band;
gamma = LFP.gamma_band;

mBetaBurst = [];mGammaEnvelope = [];
for i = 1:size(mGamma,3)
    for ii = 1:size(betaEvent{i},1)
        timestamp = floor([betaEvent{i}(ii,2)-0.075, betaEvent{i}(ii,2)+0.075].*1024);
        gammaEnvelope{i}(:,ii) = gamma(64,timestamp(1):timestamp(2));
        betaBurst{i}(:,ii) = beta(64,timestamp(1):timestamp(2));
    end
    mBetaBurst{i} = mean(betaBurst{i},2);
    mGammaEnvelope{i} = mean(gammaEnvelope{i},2);
end
mBetaBurst = horzcat(mBetaBurst{:});
mGammaEnvelope = horzcat(mGammaEnvelope{:});
%%

for i = 1:10
figure,
subplot(2,1,1),plot(mBetaBurst(:,i))
subplot(2,1,2),plot(mGammaEnvelope(:,i))
end
%%
[gammaWavelet,gf] = cwt(mean(mGammaEnvelope,2),'FrequencyLimit',[30 100],1024);
[betaWavelet,bf] = cwt(mean(mBetaBurst,2),'FrequencyLimit',[12 32],1024);
figure,
subplot(2,1,1),imagesc(-75:75,bf,abs(betaWavelet)),colormap(jet),axis xy
subplot(2,1,2),imagesc(-75:75,gf,abs(gammaWavelet)),colormap(jet),axis xy
%% Random Plots
dplot=1
if dplot
    data = mGamma(:,1);
    figure,
    subplot(311)
    plot(data)
    xlabel('Time (ms)'), ylabel('Amplitude (\muV)')
    
    % plot power (squared magnitude from origin to dot-product location in
    % complex space)
    subplot(312)
    plot(data.^2)
    xlabel('Time (ms)'), ylabel('Power \muV^2')
    
    
    % plot phase (angle of vector to dot-product, relative to positive real
    % axis)
    subplot(313)
    plot(angle(hilbert(data)))
    xlabel('Time (ms)'), ylabel('Phase (rad.)')
    
    data = mBeta(:,1);
    figure,
    subplot(311)
    plot(data)
    xlabel('Time (ms)'), ylabel('Amplitude (\muV)')
    
    % plot power (squared magnitude from origin to dot-product location in
    % complex space)
    subplot(312)
    plot(data.^2)
    xlabel('Time (ms)'), ylabel('Power \muV^2')
    
    
    % plot phase (angle of vector to dot-product, relative to positive real
    % axis)
    subplot(313)
    plot(angle(hilbert(data)))
    xlabel('Time (ms)'), ylabel('Phase (rad.)')
    
    %depth plot of gamma beta phase amplitude coupling
    figure,
    bar(betaGammaCoupling(:,30,1)),hold on
end