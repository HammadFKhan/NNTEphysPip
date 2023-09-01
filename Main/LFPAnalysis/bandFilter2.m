function LFP = bandFilter2(data,Fs,varargin)
% band pass filtering across electrodes
% also checks for 3 dimensional data and adjusts accordingly (3rd dimension
% in trials
tic
sz = size(data);
if length(sz)==3
    temp = reshape(data,sz(1),[]);
else 
    temp = data;
end
disp('Filtering band frequencies...')
init = zeros(size(temp,1),size(temp,2));
theta = init;
beta = init;
gamma = init;
parfor i = 1:size(temp,1)
    % Filter signals using FIR linear regression
    theta(i,:) = customFilt(temp(i,:),Fs,[4 10]);
    beta(i,:) = customFilt(temp(i,:),Fs,[10 30]);
    gamma(i,:) = customFilt(temp(i,:),Fs,[30 80]);
end
%% Compute Band Power
% LFP.theta_temppow  = 10*log10(abs(hilbert(theta)).^2);
% LFP.beta_temppow  = 10*log10(abs(hilbert(beta)).^2);
% LFP.gamma_temppow  = 10*log10(abs(hilbert(gamma)).^2);
LFP.theta_band = reshape(theta,sz(1),sz(2),[]);
LFP.beta_band = reshape(beta,sz(1),sz(2),[]);
LFP.gamma_band = reshape(gamma,sz(1),sz(2),[]);
toc

