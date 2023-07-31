function LFP = bandFilter(LFP,varargin)
tic
if strcmp(varargin,'depth')
    depthTrace = 1;
elseif strcmp(varargin,'single')
    depthTrace = 0;
else
    warning('defaulting filtering on best trace')
    depthTrace = 0;
end
Fs = LFP.downSampleFreq;
theta = [];
beta = [];
gamma = [];
disp('Filtering band frequencies...')

if depthTrace == 0    
    % Filter signals using FIR linear regression
    theta = customFilt(LFP.bestLFP,Fs,[4 10]);
    beta = customFilt(LFP.bestLFP,Fs,[10 30]);
    gamma = customFilt(LFP.bestLFP,Fs,[30 80]);
else
    init = zeros(size(LFP.medianLFP,1),size(LFP.medianLFP,2));
    theta = init;
    beta = init;
    gamma = init;
    temp = LFP.medianLFP;
    parfor i = 1:size(temp,1) 
        % Filter signals using FIR linear regression
        theta(i,:) = customFilt(temp(i,:),Fs,[4 10]);
        beta(i,:) = customFilt(temp(i,:),Fs,[10 30]);
        gamma(i,:) = customFilt(temp(i,:),Fs,[30 80]);
    end
end
%% Compute Band Power
LFP.theta_temppow  = 10*log10(abs(hilbert(theta)).^2);
LFP.beta_temppow  = 10*log10(abs(hilbert(beta)).^2);
LFP.gamma_temppow  = 10*log10(abs(hilbert(gamma)).^2);
LFP.theta_band = theta;
LFP.beta_band = beta;
LFP.gamma_band = gamma;
toc

