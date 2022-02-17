% Time-Frequency Analysis
function [TimeFreq,LFP,betaGroup] = tfAnalysis(Spikes,LFP)
global useGPU
tic
% Initialize Data
params.tapers = [5 9];
movingwin = [0.5 0.05];
params.pad = 2;
params.Fs = 8192;
TimeFreq = [];
% Set up Spike data
for trial = 1:size(Spikes.Clusters,2)
    fixSpiketime{trial} = Spikes.Clusters(trial).spikeTime';
end
spikeTime = sort(vertcat(fixSpiketime{:})); %% All spike times of neurons
% Set up depth specific spike data
L23 = find(Spikes.Depth.depth<250);
L4 = find(Spikes.Depth.depth>250 & Spikes.Depth.depth<400);
L5 = find(Spikes.Depth.depth>400);
% Now build spike rasters and time for each layer
clusterID = Spikes.Depth.clusterID;
l23 = {};l4 = {};l5 = {};
for i = 1:length(L23)
    l23{i} = Spikes.Clusters(clusterID(L23(i))).spikeTime'; %match electrode ID and layer location
end
for i = 1:length(L4)
    l4{i} = Spikes.Clusters(clusterID(L4(i))).spikeTime'; %match electrode ID and layer location
end
for i = 1:length(L5)
    l5{i} = Spikes.Clusters(clusterID(L5(i))).spikeTime'; %match electrode ID and layer location
end
L23spikeTime = sort(vertcat(l23{:}));
L4spikeTime = sort(vertcat(l4{:}));
L5spikeTime = sort(vertcat(l5{:}));

coherencyLFP =LFP.coherenceLFP'; % Reference LFP for spiking data (matched Fs) 
thetaLFP = LFP.theta_band'; %Downsampled LFP for P:P analysis
betaLFP = LFP.beta_band';
gammaLFP = LFP.gamma_band';

LFPTime = (0:length(coherencyLFP)-1)';
LFPTime = LFPTime/params.Fs;
downsample_LFPTime = LFP.times';

% Set up state triggered analysis
Velocity = Spikes.VR.Velocity(:,2);
velocityTrig = 4; %Triggered Velocity
loc = Spikes.VR.binWinTime*find(abs(Velocity)>velocityTrig); %multiply by the cause of bin value

% Check and adjust data structure for depth or single electrode analysis
disp('Constructing Spike Windows...')
if size(thetaLFP,1)>1
    for trial = 1:length(loc)-1
        window = [loc(trial)-1,loc(trial)+1];
        %triggered spike time with offset of the intial spike
        %this way, each spike time is translated to a window within 1 second
        timestamps(trial,:) = window; % Global timestamps for signal reference
        spike(trial).spikeTrig = spikeTime(window(1)<=spikeTime & spikeTime<=window(2))-window(1);
        % Create spike window for each layer
        L23spike(trial).spikeTrig = L23spikeTime(window(1)<=L23spikeTime & L23spikeTime<=window(2))-window(1);
        L4spike(trial).spikeTrig = L4spikeTime(window(1)<=L4spikeTime & L4spikeTime<=window(2))-window(1);
        L5spike(trial).spikeTrig = L5spikeTime(window(1)<=L5spikeTime & L5spikeTime<=window(2))-window(1);
        LFPTrig(:,:,trial) = coherencyLFP(window(1)<=LFPTime & LFPTime<=window(2));
        theta(:,:,trial)= thetaLFP(window(1)<=downsample_LFPTime & downsample_LFPTime<=window(2),:);
        beta(:,:,trial)= betaLFP(window(1)<=downsample_LFPTime & downsample_LFPTime<=window(2),:);
        gamma(:,:,trial)= gammaLFP(window(1)<=downsample_LFPTime & downsample_LFPTime<=window(2),:);
    end
else
    for trial = 1:length(loc)-1
        window = [loc(trial)-1.50,loc(trial)+.500];
        %triggered spike time with offset of the intial spike
        %this way, each spike time is translated to a window within 1 second
        timestamps(trial,:) = window; % Global timestamps for signal reference
        spike(trial).spikeTrig = spikeTime(window(1)<=spikeTime & spikeTime<=window(2))-window(1);
        LFPTrig(:,trial) = coherencyLFP(window(1)<=LFPTime & LFPTime<=window(2));
        theta(:,trial)= thetaLFP(window(1)<=downsample_LFPTime & downsample_LFPTime<=window(2));
        beta(:,trial)= betaLFP(window(1)<=downsample_LFPTime & downsample_LFPTime<=window(2));
        gamma(:,trial)= gammaLFP(window(1)<=downsample_LFPTime & downsample_LFPTime<=window(2));
    end
end

% Beta Analysis
if size(thetaLFP,2)>1
    temp = LFP.medianLFP;
    temp1 = LFP.beta_band;
    parfor electrode = 1:size(temp,1) % Checks electrode size for median
        betaGroup(electrode).electrode = groupBetaBurstDetection(temp1(electrode,:),beta(:,electrode,:),timestamps,1024); % Detect beta burst during window
    end
else
    LFP = betaBurstDetection(LFP,beta,timestamps);
end
TimeFreq.mtheta = mean(theta,2);
TimeFreq.mbeta = mean(beta,2);
LFP.betaTrials = beta;
TimeFreq.mgamma = mean(gamma,2);

%% ITPC of velocity triggered Frequency
disp('Calculating Electrode ITPC')
phaseSync = itpc(LFP,timestamps,10);
%% Oscillators Analysis using Chronux
% Phase to Phase of Theta and Beta
disp('Calculating LFP Phase Coupling')
if useGPU
    theta = gpuArray(theta);beta = gpuArray(beta);gamma = gpuArray(gamma);
params.Fs = 1024;
params.fpass = [4 30];
[tf.oscillators.tb.C,tf.oscillators.tb.phi,S12,S1,S2,tf.oscillators.t,tf.oscillators.tb.f]=cohgramc(theta,beta,movingwin,params);
params.fpass = [4 80];
[tf.oscillators.tg.C,tf.oscillators.tg.phi,S12,S1,S2,t,tf.oscillators.tg.f]=cohgramc(theta,gamma,movingwin,params);
params.fpass = [30 80];
[tf.oscillators.bg.C,tf.oscillators.bg.phi,S12,S1,S2,t,tf.oscillators.bg.f]=cohgramc(beta,gamma,movingwin,params);

params.Fs = 8192;
% Theta Band Analysis
disp('Analyzing Theta Band...')
params.fpass = [4 10];
% [tf.theta.X,tf.t,tf.theta.f]= mtspecgramc(LFPTrig,movingwin,params); 
[tf.theta.C,tf.theta.phi,S12,S1,S2,t2,tf.theta.f,zerosp]=cohgramcpt(LFPTrig,spike,movingwin,params);
% tf.mtheta = mean(tf.theta.X,3);
% tf.metheta = median(tf.theta.X,3);

params.fpass = [10 30];
% Beta Band Analysis
disp('Analyzing Beta Band...')
% [tf.beta.X,t,tf.beta.f] = mtspecgramc(LFPTrig,movingwin,params); 
[tf.beta.C,tf.beta.phi,S12,S1,S2,t2,tf.beta.f,zerosp]=cohgramcpt(LFPTrig,spike,movingwin,params); 
% tf.mbeta = mean(tf.beta.X,3);
% tf.mebeta = median(tf.beta.X,3);
params.fpass = [30 80];
% Gamma Band Analysis
disp('Analyzing Gamma Band...')
% [tf.gamma.X,t,tf.gamma.f] = mtspecgramc(LFPTrig,movingwin,params);
[tf.gamma.C,tf.gamma.phi,S12,S1,S2,t2,tf.gamma.f,zerosp]=cohgramcpt(LFPTrig,spike,movingwin,params);
% tf.mgamma = mean(tf.gamma.X,3);
% tf.megamma = median(tf.gamma.X,3);
% % Broadband
% disp('Analyzing Broad Band...')
% params.fpass = [1 80];
% [tf.broad.X,t,tf.broad.f] = mtspecgramc(LFPTrig,movingwin,params);
% [tf.broad.C,tf.broad.phi,S12,S1,S2,t2,f2,zerosp]=cohgramcpt(LFPTrig,spike,movingwin,params);
% tf.mbroad = mean(tf.broad.X,3);
% tf.mebroad = median(tf.broad.X,3);

% Depth analysis of spikes
% Theta Band Analysis
disp('Analyzing Theta Band for Layer 2/3...')
params.fpass = [4 10];
[tf.depth.L23.theta.C,tf.depth.L23.theta.phi,S12,S1,S2,t2,tf.theta.f,zerosp]=cohgramcpt(LFPTrig,L23spike,movingwin,params);
disp('Analyzing Theta Band for Layer 4...')
[tf.depth.L4.theta.C,tf.depth.L4.theta.phi,S12,S1,S2,t2,tf.theta.f,zerosp]=cohgramcpt(LFPTrig,L4spike,movingwin,params);
disp('Analyzing Theta Band for Layer 5...')
[tf.depth.L5.theta.C,tf.depth.L5.theta.phi,S12,S1,S2,t2,tf.theta.f,zerosp]=cohgramcpt(LFPTrig,L5spike,movingwin,params);

% Beta Band Analysis
params.fpass = [10 30];
disp('Analyzing Beta Band for Layer 2/3...')
[tf.depth.L23.beta.C,tf.depth.L23.beta.phi,S12,S1,S2,t2,tf.beta.f,zerosp]=cohgramcpt(LFPTrig,L23spike,movingwin,params);
disp('Analyzing Beta Band for Layer 4...')
[tf.depth.L4.beta.C,tf.depth.L4.beta.phi,S12,S1,S2,t2,tf.beta.f,zerosp]=cohgramcpt(LFPTrig,L4spike,movingwin,params);
disp('Analyzing Beta Band for Layer 5...')
[tf.depth.L5.beta.C,tf.depth.L5.beta.phi,S12,S1,S2,t2,tf.beta.f,zerosp]=cohgramcpt(LFPTrig,L5spike,movingwin,params);

% Gamma Band Analysis
params.fpass = [30 80];
disp('Analyzing Gamma Band for Layer 2/3...')
[tf.depth.L23.gamma.C,tf.depth.L23.gamma.phi,S12,S1,S2,t2,tf.gamma.f,zerosp]=cohgramcpt(LFPTrig,L23spike,movingwin,params);
disp('Analyzing Gamma Band for Layer 4...')
[tf.depth.L4.gamma.C,tf.depth.L4.gamma.phi,S12,S1,S2,t2,tf.gamma.f,zerosp]=cohgramcpt(LFPTrig,L4spike,movingwin,params);
disp('Analyzing Gamma Band for Layer 5...')
[tf.depth.L5.gamma.C,tf.depth.L5.gamma.phi,S12,S1,S2,t2,tf.gamma.f,zerosp]=cohgramcpt(LFPTrig,L5spike,movingwin,params);

%% Inter-trial Phase Clustering of Spike-LFP
[tf.theta.theta,tf.theta.itpc,tf.beta.beta,tf.beta.itpc,tf.gamma.gamma,tf.gamma.itpc] = chronuxITPC(tf.theta.phi,tf.beta.phi,tf.gamma.phi);

%% Inter-trial Phase Clustering of Spike-LFP per Layer
disp('ITPC for layer 2/3')
[tf.depth.L23.theta.theta,tf.depth.L23.theta.itpc,tf.depth.L23.beta.beta,tf.depth.L23.beta.itpc,tf.depth.L23.gamma.gamma,tf.depth.L23.gamma.itpc]...
    = chronuxITPC(tf.depth.L23.theta.phi,tf.depth.L23.beta.phi,tf.depth.L23.gamma.phi);
disp('ITPC for layer 4')
[tf.depth.L4.theta.theta,tf.depth.L4.theta.itpc,tf.depth.L4.beta.beta,tf.depth.L4.beta.itpc,tf.depth.L4.gamma.gamma,tf.depth.L4.gamma.itpc]...
    = chronuxITPC(tf.depth.L4.theta.phi,tf.depth.L4.beta.phi,tf.depth.L4.gamma.phi);
disp('ITPC for layer 5')
[tf.depth.L5.theta.theta,tf.depth.L5.theta.itpc,tf.depth.L5.beta.beta,tf.depth.L5.beta.itpc,tf.depth.L5.gamma.gamma,tf.depth.L5.gamma.itpc]...
    = chronuxITPC(tf.depth.L5.theta.phi,tf.depth.L5.beta.phi,tf.depth.L5.gamma.phi);


% disp('Performing Wavelet Analysis...')
% for trial = 1:size(LFPTrig,2)
%     [tf.wavelet.cfs(:,:,trial),tf.wavelet.f] = cwt(LFPTrig(:,trial),1024,'FrequencyLimits',[0 80]);
%     [tf.wavelet.theta_cfs(:,:,trial),tf.wavelet.theta_f] = cwt(theta(:,trial),1024,'FrequencyLimits',[4 10]);
%     
%     [tf.wavelet.beta_cfs(:,:,trial),tf.wavelet.beta_f] = cwt(beta(:,trial),1024,'FrequencyLimits',[10 30]);
%     
%     [tf.wavelet.gamma_cfs(:,:,trial),tf.wavelet.gamma_f] = cwt(gamma(:,trial),1024,'FrequencyLimits',[30 80]);
% 
% end
TimeFreq.tf = tf;
toc
end
