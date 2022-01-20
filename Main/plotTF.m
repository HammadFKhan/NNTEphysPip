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

% Spike-coherence
figure,imagesc(-1000:1000,tf.theta.f,interp2(mean(tf.theta.C,3),interpLevel)),caxis([0 .5]),title('\theta Spike Coupling'),colormap(jet),axis xy,set(gcf, 'Position',  [100, 100, 500, 500])
figure,imagesc(-1000:1000,tf.beta.f,interp2(mean(tf.beta.C,3),interpLevel)),caxis([0 .5]),title('\beta Spike Coupling'),colormap(jet),axis xy,set(gcf, 'Position',  [100, 100, 500, 500])
figure,imagesc(-1000:1000,tf.gamma.f,interp2(mean(tf.gamma.C,3),interpLevel)),caxis([0 .5]),title('\gamma Spike Coupling'),colormap(jet),axis xy,set(gcf, 'Position',  [100, 100, 500, 500])


figure('Name','All Layers'),
% plotpad = size(tf.theta.theta,2);
% polarplot([zeros(1,plotpad); tf.theta.theta],[zeros(1,plotpad); ones(1,plotpad)],'k');
subplot(3,1,1),polarhistogram(horzcat(tf.theta.theta{:}));
title(['Theta ITPC: ' num2str(mean(tf.theta.itpc))]);
subplot(3,1,2),polarhistogram(horzcat(tf.beta.beta{:}),20);
title(['Beta ITPC: ' num2str(mean(tf.beta.itpc))]);
subplot(3,1,3),polarhistogram(horzcat(tf.gamma.gamma{:}),20);
title(['Gamma ITPC: ' num2str(mean(tf.gamma.itpc))]);
% Plots by depth
if isfield(TimeFreq.tf,'depth')
    
    % Plot Spike-coherence
    % L2/3
    % check for nan values
    keep = [];
    check = [];
    check = isnan(tf.depth.L23.theta.C);
    [~,~,t] = findND(check==1); % Find trials that have nans
    t = unique(t);
    keep = 1:size(tf.depth.L23.theta.C,3);
    keep(:,t) = [];
    
    figure,imagesc(-1000:1000,tf.theta.f,interp2(mean(tf.depth.L23.theta.C(:,:,keep),3),interpLevel))...
        ,caxis([0 .5]),title('\theta-L2/3 Spike Coupling'),colormap(jet),axis xy,set(gcf, 'Position',  [100, 100, 500, 500])
    figure,imagesc(-1000:1000,tf.beta.f,interp2(median(tf.depth.L23.beta.C(:,:,keep),3),interpLevel))...
        ,caxis([0 .5]),title('\beta-L2/3 Spike Coupling'),colormap(jet),axis xy,set(gcf, 'Position',  [100, 100, 500, 500])
    figure,imagesc(-1000:1000,tf.gamma.f,interp2(median(tf.depth.L23.gamma.C(:,:,keep),3),interpLevel))...
        ,caxis([0 .5]),title('\gamma-L2/3 Spike Coupling'),colormap(jet),axis xy,set(gcf, 'Position',  [100, 100, 500, 500])
    
    % L4
    % check for nan values
    keep = [];
    check = [];
    check = isnan(tf.depth.L4.theta.C);
    [~,~,t] = findND(check==1); % Find trials that have nans
    t = unique(t);
    keep = 1:size(tf.depth.L4.theta.C,3);
    keep(:,t) = [];
    figure,imagesc(-1000:1000,tf.theta.f,interp2(median(tf.depth.L4.theta.C(:,:,keep),3),interpLevel))...
        ,caxis([0 .5]),title('\theta-L4 Spike Coupling'),colormap(jet),axis xy,set(gcf, 'Position',  [100, 100, 500, 500])
    figure,imagesc(-1000:1000,tf.beta.f,interp2(median(tf.depth.L4.beta.C(:,:,keep),3),interpLevel))...
        ,caxis([0 .5]),title('\beta-L4 Spike Coupling'),colormap(jet),axis xy,set(gcf, 'Position',  [100, 100, 500, 500])
    figure,imagesc(-1000:1000,tf.gamma.f,interp2(median(tf.depth.L4.gamma.C(:,:,keep),3),interpLevel))...
        ,caxis([0 .5]),title('\gamma-L4 Spike Coupling'),colormap(jet),axis xy,set(gcf, 'Position',  [100, 100, 500, 500])
    % L5
    % check for nan values
    keep = [];
    check = [];
    check = isnan(tf.depth.L5.theta.C);
    [~,~,t] = findND(check==1); % Find trials that have nans
    t = unique(t);
    keep = 1:size(tf.depth.L5.theta.C,3);
    keep(:,t) = [];
    
    figure,imagesc(-1000:1000,tf.theta.f,interp2(mean(tf.depth.L5.theta.C(:,:,keep),3),interpLevel))...
        ,caxis([0 .5]),title('\theta-L5 Spike Coupling'),colormap(jet),axis xy,set(gcf, 'Position',  [100, 100, 500, 500])
    figure,imagesc(-1000:1000,tf.beta.f,interp2(mean(tf.depth.L5.beta.C(:,:,keep),3),interpLevel))...
        ,caxis([0 .5]),title('\beta-L5 Spike Coupling'),colormap(jet),axis xy,set(gcf, 'Position',  [100, 100, 500, 500])
    figure,imagesc(-1000:1000,tf.gamma.f,interp2(mean(tf.depth.L5.gamma.C(:,:,keep),3),interpLevel))...
        ,caxis([0 .5]),title('\gamma-L5 Spike Coupling'),colormap(jet),axis xy,set(gcf, 'Position',  [100, 100, 500, 500])
    
    
    % Plot ITPC
    figure('Name','Layer2/3'),
    subplot(3,1,1),polarhistogram(horzcat(tf.depth.L23.theta.theta{:}));
    title(['Theta ITPC: ' num2str(mean(tf.depth.L23.theta.itpc))]);
    subplot(3,1,2),polarhistogram(horzcat(tf.depth.L23.beta.beta{:}),20);
    title(['Beta ITPC: ' num2str(mean(tf.depth.L23.beta.itpc))]);
    subplot(3,1,3),polarhistogram(horzcat(tf.depth.L23.gamma.gamma{:}),20);
    title(['Gamma ITPC: ' num2str(mean(tf.depth.L23.gamma.itpc))]);
    
    figure('Name','Layer4'),
    subplot(3,1,1),polarhistogram(horzcat(tf.depth.L4.theta.theta{:}));
    title(['Theta ITPC: ' num2str(mean(tf.depth.L4.theta.itpc))]);
    subplot(3,1,2),polarhistogram(horzcat(tf.depth.L4.beta.beta{:}),20);
    title(['Beta ITPC: ' num2str(mean(tf.depth.L4.beta.itpc))]);
    subplot(3,1,3),polarhistogram(horzcat(tf.depth.L4.gamma.gamma{:}),20);
    title(['Gamma ITPC: ' num2str(mean(tf.depth.L4.gamma.itpc))]);
    
    figure('Name','Layer5'),
    subplot(3,1,1),polarhistogram(horzcat(tf.depth.L5.theta.theta{:}));
    title(['Theta ITPC: ' num2str(mean(tf.depth.L5.theta.itpc))]);
    subplot(3,1,2),polarhistogram(horzcat(tf.depth.L5.beta.beta{:}),20);
    title(['Beta ITPC: ' num2str(mean(tf.depth.L5.beta.itpc))]);
    subplot(3,1,3),polarhistogram(horzcat(tf.depth.L5.gamma.gamma{:}),20);
    title(['Gamma ITPC: ' num2str(mean(tf.depth.L5.gamma.itpc))]);
end

figure('Name','All Layers'),
itpcAll = vertcat(tf.theta.itpc,tf.beta.itpc,tf.gamma.itpc);
itpcF = vertcat(tf.theta.f',tf.beta.f',tf.gamma.f');
plot(itpcF,itpcAll);
% For each layer
figure('Name','L2/3')
itpcL23 = vertcat(tf.depth.L23.theta.itpc,tf.depth.L23.beta.itpc,tf.depth.L23.gamma.itpc);
plot(itpcF,itpcL23);
figure('Name','L4')
itpcL4 = vertcat(tf.depth.L4.theta.itpc,tf.depth.L4.beta.itpc,tf.depth.L4.gamma.itpc);
plot(itpcF,itpcL4);
figure('Name','L5')
itpcL5 = vertcat(tf.depth.L5.theta.itpc,tf.depth.L5.beta.itpc,tf.depth.L5.gamma.itpc);
plot(itpcF,itpcL5);
end