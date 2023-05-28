function output = gammaBetaCoupling(LFP,Tf,betaGroup)
% Analyze coupling stats between beta/gamma LFPs and beta/gamma events

mGamma = Tf.gammaLFP;
mBeta = Tf.betaLFP;
%% Gamma amplitude phase coupling by depth and trial
fprintf('Calculating beta-gamma amplitude coupling...')
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
%         bar(-179:20:180,betaGammaCoupling)
    end
%     xlim([-180 180])
end
fprintf('done\n')
%% gamma beta LFP coupling during behavior state
betaLFP = Tf.betaLFP;
gammaLFP = Tf.gammaLFP;
betaLFP = permute(betaLFP,[2 3 1]);
gammaLFP = permute(gammaLFP,[2 3 1]);

betaDepth = mean(betaLFP.^2,3);
gammaDepth = mean(gammaLFP.^2,3);
betaStd = std(betaDepth');
gammaStd = std(gammaDepth');

betaDepth = mean(betaDepth,2);
gammaDepth = mean(gammaDepth,2);

betaDepthAd = betaDepth-min(betaDepth,[],'all');
gammaDepthAd = (betaDepth-gammaDepth);gammaDepthAd = gammaDepthAd-min(gammaDepthAd,[],'all');
buff = betaStd-gammaStd; %betaStd-gammaStd;
betaDepthNorm = (betaDepthAd-min(betaDepthAd,[],'all'))/(max(betaDepthAd,[],'all')-min(betaDepthAd,[],'all'));
gammaDepthNorm = (gammaDepthAd-min(gammaDepthAd,[],'all'))/(max(gammaDepthAd,[],'all')-min(gammaDepthAd,[],'all'));
betaStdNorm = (betaStd-min(betaStd,[],'all'))/(max(betaStd,[],'all')-min(betaStd,[],'all'));
gammaStdNorm = (buff-min(buff,[],'all'))/(max(buff,[],'all')-min(buff,[],'all'));


%% beta event/gamma amplitude modulation
% gamma amplitude during beta events
mBetaBurst = [];mGammaWaveform = [];mGammaEnvelope = [];mBetaPower = [];
for electrode = 1:64
    fprintf('Analyzing electrode: %d...',electrode')
    betaEvent = betaGroup(electrode).electrode.betaBurst.detectedBeta; %use global LFP values as beta bursts is recorded with this reference
    beta = LFP.beta_band;
    gamma = LFP.gamma_band;
    for i = 1:size(mGamma,3)
        for ii = 1:size(betaEvent{i},1)
            timestamp = floor([betaEvent{i}(ii,2)-0.075, betaEvent{i}(ii,2)+0.075].*1024);
            gammaWaveform{i}(:,ii) = gamma(electrode,timestamp(1):timestamp(2));
            betaBurst{i}(:,ii) = beta(electrode,timestamp(1):timestamp(2));
        end
        if isempty(ii)
            BetaBurst{i} = [];
            GammaWaveform{i} = [];
        else
            BetaBurst{i} = mean(betaBurst{i},2);
            BetaPower{i} = mean(betaBurst{i}.^2,2);
            GammaWaveform{i} = mean(gammaWaveform{i},2);
            GammaEnvelope{i} = mean(gammaWaveform{i}.^2,2);
        end
    end
    mBetaBurst{electrode} = horzcat(BetaBurst{:}); % Beta burst
    mBetaPower{electrode} = horzcat(BetaPower{:}); % Beta power
    mGammaWaveform{electrode} = horzcat(GammaWaveform{:}); %Gamma LFP
    mGammaEnvelope{electrode} = horzcat(GammaEnvelope{:}); %Gamma amplitude
    fprintf('done\n')
end
% beta event/gamme amplitude distribution across electrodes
% mBetaBurstDepth = cellfun(@(mBetaBurst) mean(mBetaBurst,2),mBetaBurst,'UniformOutput',false);
% mBetaPowerDepth = cellfun(@(mBetaPower) mean(mBetaPower,2),mBetaPower,'UniformOutput',false);
% mGammaWaveformDepth = cellfun(@(mGammaWaveform) mean(mGammaWaveform,2), mGammaWaveform, 'UniformOutput',false);
% mGammaEnvelopeDepth = cellfun(@(mGammaEnvelope) mean(mGammaEnvelope,2), mGammaEnvelope, 'UniformOutput', false);
mBetaBurstDepth = cellfun(@mean,mBetaBurst,'UniformOutput',false);
mBetaPowerDepth = cellfun(@mean,mBetaPower,'UniformOutput',false);
mGammaWaveformDepth = cellfun(@mean,mGammaWaveform, 'UniformOutput',false);
mGammaEnvelopeDepth = cellfun(@mean,mGammaEnvelope, 'UniformOutput', false);
%calculate mean and std
betaPowerDepth = cellfun(@mean,mBetaPowerDepth);
betaPowerStd = cellfun(@std,mBetaPowerDepth);
gammaEnvelopeDepth = cellfun(@mean,mGammaEnvelopeDepth);
gammaEnvelopeStd = cellfun(@std,mGammaEnvelopeDepth);

%Normalize data
betaEventDepthAd = betaPowerDepth-min(betaPowerDepth,[],'all');
gammaEnvelopeDepthAd = (betaPowerDepth-gammaEnvelopeDepth);
gammaEnvelopeDepthAd = gammaEnvelopeDepthAd-min(gammaEnvelopeDepthAd,[],'all');
buff = gammaEnvelopeStd; %betaStd-gammaStd;
betaEventDepthNorm = (betaEventDepthAd-min(betaEventDepthAd,[],'all'))...
    /(max(betaEventDepthAd,[],'all')-min(betaEventDepthAd,[],'all'));
gammaEnvelopeDepthNorm = (gammaEnvelopeDepthAd-min(gammaEnvelopeDepthAd,[],'all'))...
    /(max(gammaEnvelopeDepthAd,[],'all')-min(gammaEnvelopeDepthAd,[],'all'));
betaEventStdNorm = (betaPowerStd-min(betaStd,[],'all'))/(max(betaPowerStd,[],'all')-min(betaPowerStd,[],'all'));
gammaEnvelopeStdNorm = (buff-min(buff,[],'all'))/(max(buff,[],'all')-min(buff,[],'all'));
%% Output
output.betaDepthNorm = betaDepthNorm; 
output.betaDepthStd = betaStdNorm;
output.gammaDepthNorm = gammaDepthNorm;
output.gammaDepthStd = gammaStdNorm;
output.mBetaBurst = mBetaBurst;

output.mBetaPowerDepth = mBetaPowerDepth;
output.mGammaWaveform = mGammaWaveform;
output.mGammaEnvelope = mGammaEnvelope;

output.betaGammaCoupling = betaGammaCoupling;
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
    
    % gamma envelope beta event plot by depth
%     plotBeta = cellfun(@(mBetaPower) mean(mBetaPower,2),mBetaPower,'UniformOutput',false); % format for plotting 
%     plotGamma = cellfun(@(mGammaEnvelope) mean(mGammaEnvelope,2),mGammaEnvelope,'UniformOutput',false);
%     figure,
%     subplot(2,1,1),stack_plot(horzcat(plotBeta{:})',1,5);
%     subplot(2,1,2),stack_plot(horzcat(plotGamma{:})',1,25);
    
    % wavelet decomposition of plots
    plotBeta = cellfun(@(mBetaBurst) mean(mBetaBurst,2),mBetaBurst,'UniformOutput',false); % format for plotting 
    plotGamma = cellfun(@(mGammaWaveform) mean(mGammaWaveform,2),mGammaWaveform,'UniformOutput',false);
    plotBeta = mean(horzcat(plotBeta{:}),2);
    plotGamma = mean(horzcat(plotGamma{:}),2);
    [gammaWavelet,gf] = cwt(plotBeta,'FrequencyLimit',[30 100],1024);
    [betaWavelet,bf] = cwt(plotBeta,'FrequencyLimit',[12 32],1024);
    figure,
    subplot(2,1,1),imagesc(-75:75,bf,abs(betaWavelet)),colormap(jet),axis xy,hold on
    subplot(2,1,2),imagesc(-75:75,gf,abs(gammaWavelet)),colormap(jet),axis xy
    figure,
    subplot(2,1,1),plot(plotBeta)
    subplot(2,1,2),plot(plotGamma)
    % Gamma beta LFP profile plots
    betaData = betaDepthNorm; betaDataStd= betaStdNorm'; gammaData = gammaDepthNorm; gammaDataStd = gammaStdNorm';
    figure,
    subplot(2,1,1),plot(1:64,smoothdata(betaData,'gaussian',10)),hold on,title('Beta LFP')
    subplot(2,1,1),plot(1:64,smoothdata(betaData+betaDataStd,'gaussian',10),'r--');
    subplot(2,1,1),plot(1:64,smoothdata(betaData-betaDataStd,'gaussian',10),'r--'),ylim([0. .9]),box off,set(gca,'TickDir','out');
    
    subplot(2,1,2),plot(1:64,smoothdata(gammaData,'gaussian',10)),hold on,title('Gamma LFP')
    subplot(2,1,2),plot(1:64,smoothdata(gammaData+gammaDataStd,'gaussian',10),'r--');
    subplot(2,1,2),plot(1:64,smoothdata(gammaData-gammaDataStd,'gaussian',10),'r--'),ylim([0. .9]),box off,set(gca,'TickDir','out');
    
    % Gamma beta Event profile
%     figure, %non adjusted data
%     subplot(2,1,1),plot(smoothdata(betaPowerDepth,'gaussian',10))
%     subplot(2,1,2),plot(smoothdata(gammaEnvelopeDepth,'gaussian',10))
    
    figure,
    subplot(2,1,1),plot(1:64,smoothdata(betaEventDepthNorm,'gaussian',10)),hold on,title('Beta Event')
    subplot(2,1,1),plot(1:64,smoothdata(betaEventDepthNorm+betaEventStdNorm,'gaussian',10),'r--');
    subplot(2,1,1),plot(1:64,smoothdata(betaEventDepthNorm-betaEventStdNorm,'gaussian',10),'r--'),ylim([0. .9]),box off,set(gca,'TickDir','out');
    
    subplot(2,1,2),plot(1:64,smoothdata(gammaEnvelopeDepthNorm,'gaussian',10)),hold on,title('Gamma Envelope')
    subplot(2,1,2),plot(1:64,smoothdata(gammaEnvelopeDepthNorm+gammaEnvelopeStdNorm,'gaussian',10),'r--');
    subplot(2,1,2),plot(1:64,smoothdata(gammaEnvelopeDepthNorm-gammaEnvelopeStdNorm,'gaussian',10),'r--'),ylim([0. .9]),box off,set(gca,'TickDir','out');

end