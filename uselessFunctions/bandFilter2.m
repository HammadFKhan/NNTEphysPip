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
    theta(i,:) = customFilt(temp(i,:),Fs,[4 12]);
    beta(i,:) = customFilt(temp(i,:),Fs,[12 30]);
    gamma(i,:) = customFilt(temp(i,:),Fs,[30 90]);
end
% for i = 1:size(temp,1)
%     thetaPow(i,:) = abs(hilbert(theta(i,:)));
%     betaPow(i,:) = abs(hilbert(beta(i,:)));
%     gammaPow(i,:) = abs(hilbert(gamma(i,:)));
% end
% %% Compute Band Power
% LFP.thetaPower  = reshape(thetaPow,sz(1),sz(2),[]);
% LFP.betaPower  = reshape(betaPow,sz(1),sz(2),[]);
% LFP.gammaPower = reshape(gammaPow,sz(1),sz(2),[]);
LFP.theta_band = reshape(theta,sz(1),sz(2),[]);
LFP.beta_band = reshape(beta,sz(1),sz(2),[]);
LFP.gamma_band = reshape(gamma,sz(1),sz(2),[]);
toc

