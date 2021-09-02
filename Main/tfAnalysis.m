% Time-Frequency Analysis
function TimeFreq = tfAnalysis(Spikes,LFP)
% Initialize Data
params.tapers = [5 9];
movingwin = [0.5 0.05];
params.pad = 2;
params.Fs = 8192;

% Set up Spike data
for trial = 1:size(Spikes.Clusters,2)
    fixSpiketime{trial} = Spikes.Clusters(trial).spikeTime';
end
spikeTime = sort(vertcat(fixSpiketime{:}));

coherencyLFP =LFP.coherenceLFP'; % Reference LFP for spiking data
thetaLFP = LFP.theta_band'; %Downsampled LFP for P:P analysis
betaLFP = LFP.beta_band';
gammaLFP = LFP.gamma_band';

LFPTime = (0:length(coherencyLFP)-1)';
LFPTime = LFPTime/params.Fs;
downsample_LFPTime = LFP.times';

% Set up state triggered analysis
Velocity = Spikes.VR.Velocity(:,2);
velocityTrig = 1.2*mean(Velocity);
loc = 3*find(Velocity>velocityTrig); %multiply by the cause of bin value
for trial = 1:length(loc)-1
    window = (loc(trial)-1.00:loc(trial)+1.00);
    %triggered spike time with offset of the intial spike
    %this way, each spike time is translated to a window within 1 second
    spike(trial).spikeTrig = spikeTime(window(1)<=spikeTime & spikeTime<=window(2))-window(1);
    LFPTrig(:,trial) = coherencyLFP(window(1)<=LFPTime & LFPTime<=window(2));
    theta(:,trial)= thetaLFP(window(1)<=downsample_LFPTime & downsample_LFPTime<=window(2));
    beta(:,trial)= betaLFP(window(1)<=downsample_LFPTime & downsample_LFPTime<=window(2));
    gamma(:,trial)= gammaLFP(window(1)<=downsample_LFPTime & downsample_LFPTime<=window(2));

end
TimeFreq.mtheta = mean(theta,2);
TimeFreq.mbeta = mean(beta,2);
TimeFreq.mgamma = mean(gamma,2);
% Oscillators Analysis
% Phase to Phase of Theta and Beta
params.Fs = 1024;
params.fpass = [4 30];
[tf.oscillators.tb.C,tf.oscillators.tb.phi,S12,S1,S2,tf.oscillators.t,tf.oscillators.tb.f]=cohgramc(theta,beta,movingwin,params);
params.fpass = [4 80];
[tf.oscillators.tg.C,tf.oscillators.tg.phi,S12,S1,S2,t,tf.oscillators.tg.f]=cohgramc(theta,gamma,movingwin,params);
params.fpass = [30 80];
[tf.oscillators.bg.C,tf.oscillators.bg.phi,S12,S1,S2,t,tf.oscillators.bg.f]=cohgramc(beta,gamma,movingwin,params);

params.Fs = 8192;
% Theta Band Analysis
% disp('Analyzing Theta Band...')
% params.fpass = [4 10];
% [tf.theta.X,tf.t,tf.theta.f]= mtspecgramc(LFPTrig,movingwin,params); 
% [tf.theta.C,tf.theta.phi,S12,S1,S2,t2,f2,zerosp]=cohgramcpt(LFPTrig,spike,movingwin,params);
% tf.mtheta = mean(tf.theta.X,3);
% tf.metheta = median(tf.theta.X,3);

% params.fpass = [10 30];
% % Beta Band Analysis
% disp('Analyzing Beta Band...')
% [tf.beta.X,t,tf.beta.f] = mtspecgramc(LFPTrig,movingwin,params); 
% [tf.beta.C,tf.beta.phi,S12,S1,S2,t2,f2,zerosp]=cohgramcpt(LFPTrig,spike,movingwin,params); 
% tf.mbeta = mean(tf.beta.X,3);
% tf.mebeta = median(tf.beta.X,3);
% params.fpass = [30 80];
% % Gamma Band Analysis
% disp('Analyzing Gamma Band...')
% [tf.gamma.X,t,tf.gamma.f] = mtspecgramc(LFPTrig,movingwin,params);
% [tf.gamma.C,tf.gamma.phi,S12,S1,S2,t2,f2,zerosp]=cohgramcpt(LFPTrig,spike,movingwin,params);
% tf.mgamma = mean(tf.gamma.X,3);
% tf.megamma = median(tf.gamma.X,3);
% % Broadband
% disp('Analyzing Broad Band...')
% params.fpass = [1 80];
% [tf.broad.X,t,tf.broad.f] = mtspecgramc(LFPTrig,movingwin,params);
% [tf.broad.C,tf.broad.phi,S12,S1,S2,t2,f2,zerosp]=cohgramcpt(LFPTrig,spike,movingwin,params);
% tf.mbroad = mean(tf.broad.X,3);
% tf.mebroad = median(tf.broad.X,3);
% Inter-trial Phase Clustering
% tf.theta.theta = tf.theta.phi(:,2,:);
% tf.beta.beta = tf.beta.phi(:,2,:);
% tf.gamma.gamma = tf.gamma.phi(:,2,:);
% 
% tf.theta.itpc = abs(mean(exp(1i*tf.theta.theta)));
% tf.beta.itpc = abs(mean(exp(1i*tf.beta.beta)));
% tf.gamma.itpc = abs(mean(exp(1i*tf.gamma.gamma)));
disp('Performing Wavelet Analysis...')
for trial = 1:size(LFPTrig,2)
    [tf.wavelet.cfs(:,:,trial),tf.wavelet.f] = cwt(LFPTrig(:,trial),1024,'FrequencyLimits',[0 80]);
    [tf.wavelet.theta_cfs(:,:,trial),tf.wavelet.theta_f] = cwt(theta(:,trial),1024,'FrequencyLimits',[4 10]);
    
    [tf.wavelet.beta_cfs(:,:,trial),tf.wavelet.beta_f] = cwt(beta(:,trial),1024,'FrequencyLimits',[10 30]);
    
    [tf.wavelet.gamma_cfs(:,:,trial),tf.wavelet.gamma_f] = cwt(gamma(:,trial),1024,'FrequencyLimits',[30 80]);

end
TimeFreq.tf = tf;
end
