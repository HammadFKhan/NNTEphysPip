%% Ephys from Austin
if ~exist('trialMatrix','var')
   % load('Y:\Austin\ActiveTouchExperiment\240614_M2ActiveTouch\Day2\TrialLFPData.mat')
    load('Y:\Austin\ActiveTouchExperiment\240614_M2ActiveTouch\Day1\TrialLFPData.mat')
end
trialMatrix = permute(trialMatrix,[2 3 1]);
%%
sz = size(trialMatrix);
xo = [];
for n = 1:sz(3)
    dat = trialMatrix(:,:,n);
    dat(dat==0) = dat(dat==0)+0.0001;
    dat(isnan(dat)) = 0;
     a = squeeze(mean(dat,1));
    [PSDCh(n,:),f] = pwelch(a,300,0,1000,FsG);
    xo(:,:,n) = bandpass_filter(dat,5,40,FsG); %x,f1,f2,Fs
end
%%
time = 1/FsG:1/FsG:sz(2)/FsG;
figure
plot(time(1:10:end),(mean(xo(:,1:10:end,:),3)'),'color',[0.5 0.5 0.5]);
%%
figure,
semilogy(f,mean(PSDCh)),xlim([0 100])
%% 
if ~exist('trialStuct','var')
   % load('Y:\Austin\ActiveTouchExperiment\240614_M2ActiveTouch\Day2\TrialLFPData.mat')
    load('Y:\Austin\ActiveTouchExperiment\240614_M2ActiveTouch\Day1\TrialLFPData.mat')
end
t= arrayfun(@(x) cell2mat(x.IndividualTrial), trialStruct, 'UniformOutput', false);
for n = 1:162
    trialStruct(n).IndividualTrial = t{n};
end
%%
