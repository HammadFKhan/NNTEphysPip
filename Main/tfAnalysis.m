% Time-Frequency Analysis
function Spikes = tfAnalysis(Spikes,LFP)
% Initialize Data
params.Fs = 8192;
params.tapers = [5 9];
movingwin = [0.5 0.05];
params.pad = 2;
% Set up Spike data
spikeData = Spikes.Clusters(1).cluster;
data = LFP';
for i = 1:size(Spikes.Clusters,2)
    fixSpiketime{i} = Spikes.Clusters(i).spikeTime';
end
spikeTime = sort(vertcat(fixSpiketime{:}));

coherencyLFP =LFP'; % Reference LFP (need to change later)
LFPTime = (1:length(coherencyLFP))';
LFPTime = LFPTime/params.Fs;

% Set up state triggered analysis
Velocity = Spikes.VR.Velocity(:,2);
velocityTrig = 1.2*mean(Velocity);
loc = 3*find(Velocity>velocityTrig); %multiply by the cause of bin value
for i = 1:length(loc)-1
    window = (loc(i)-0.500:loc(i)+.500);
    spike(i).spikeTrig = spikeTime(window(1)<=spikeTime & spikeTime<=window(2));
    LFPTrig(:,i) = coherencyLFP(window(1)<=LFPTime & LFPTime<=window(2));
end
% Theta Band Analysis
disp('Analyzing Theta Band...')
params.fpass = [0 10];
[tf.theta.X,tf.t,tf.theta.f]= mtspecgramc(LFPTrig,movingwin,params); 
[tf.theta.C,tf.theta.phi,S12,S1,S2,t2,f2,zerosp]=cohgramcpt(LFPTrig,spike,movingwin,params);
tf.mtheta = mean(tf.theta.X,3);
tf.metheta = median(tf.theta.X,3);
params.fpass = [10 30];
% Beta Band Analysis
disp('Analyzing Beta Band...')
[tf.beta.X,t,tf.beta.f] = mtspecgramc(data,movingwin,params); 
[tf.beta.C,tf.beta.phi,S12,S1,S2,t2,f2,zerosp]=cohgramcpt(LFPTrig,spike,movingwin,params); 
tf.mbeta = mean(tf.beta.X,3);
tf.mebeta = median(tf.beta.X,3);
params.fpass = [30 80];
% Gamma Band Analysis
disp('Analyzing Gamma Band...')
[tf.gamma.X,t,tf.gamma.f] = mtspecgramc(data,movingwin,params);
[tf.gamma.C,tf.gamma.phi,S12,S1,S2,t2,f2,zerosp]=cohgramcpt(LFPTrig,spike,movingwin,params);
tf.mgamma = mean(tf.gamma.X,3);
tf.megamma = median(tf.gamma.X,3);
% Broadband
% disp('Analyzing Broad Band...')
% params.fpass = [0 80];
% [tf.broad.X,t,tf.broad.f] = mtspecgramc(data,movingwin,params);
% [tf.broad.C,tf.broad.phi,S12,S1,S2,t2,f2,zerosp]=cohgramcpt(coherencyLFP,spikeTime,movingwin,params);

% Inter-trial Phase Clustering
% tf.theta.theta = tf.theta.phi(:,2);
% tf.beta.beta = tf.beta.phi(:,2);
% tf.gamma.gamma = tf.gamma.phi(:,2);
% 
% tf.theta.itpc = abs(mean(exp(1i*tf.theta.theta)));
% tf.beta.itpc = abs(mean(exp(1i*tf.beta.beta)));
% tf.gamma.itpc = abs(mean(exp(1i*tf.gamma.gamma)));

Spikes.tf = tf;
end
