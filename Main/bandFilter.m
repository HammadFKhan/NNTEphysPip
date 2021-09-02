function LFP = bandFilter(LFP)

Fs = LFP.downSampleFreq;
% Filter signals using FIR linear regression
disp('Filtering band frequencies...')
theta = customFilt(LFP.bestLFP,Fs,[4 10]);
beta = customFilt(LFP.bestLFP,Fs,[10 30]);
gamma = customFilt(LFP.bestLFP,Fs,[30 80]);
%% Compute Band Power
LFP.theta_temppow  = 10*log10(abs(hilbert(theta)).^2);
LFP.beta_temppow  = 10*log10(abs(hilbert(beta)).^2);
LFP.gamma_temppow  = 10*log10(abs(hilbert(gamma)).^2);

LFP.theta_band = theta;
LFP.beta_band = beta;
LFP.gamma_band = gamma;
