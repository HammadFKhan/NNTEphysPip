%% Wave data analysis 
% BdChnl is the grid index (before remapping) that has a resistance >10mOhm
% ChanSpace is the grid spacing, fb, ff, and sw are the wave properties
% (start, source, speed, length, duration, direction, and evaluation
% point). Each cell is one of the properties I don’t remember the order I’d
% have to check. ff contains waves >0 and <100ms post stimulus, fb is >100
% and <300 post stim, sw is all other waves. StimData is the lfp 1.5s
% before and after stimulus, the dimensions are trails x channel x time.
% Thr is the threshold for determining if something is a wave, which is
% determined by shuffling the channels and calculating rho, then finding
% the 99% of rhos
clear
[file,fpath] = uigetfile('Z:\data\Austin\GridData');
load(fullfile(fpath,file))
%%% Remake data to structures for readability
Waves.ff = remakeDat(ff);
Waves.fb = remakeDat(fb);
Waves.sw = remakeDat(sw);
%%% save data in fpath
sessionName = [fpath,'/','Waves.mat'];
fprintf('Saving data...')
save(sessionName,"Waves","BdChnl","fpath","stimData","thr","ChanSpace","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
fprintf('done!\n')
%% Analyze wave properties and LFP 
% TODO: plot Phase alignment 
if ~exist('Waves','var')
load('Z:\data\Austin\GridData\240508_eOPN_Inhib\WavesCombined.mat')
%load('Z:\data\Austin\GridData\240523_eOPNInhib\WavesCombined.mat')
end


%% FEED FORWARD
f = figure(1);
dat1 = Waves.baseline.ff.speed;
dat2 = Waves.eOPN.ff.speed;
temp = nan(max([length(dat1) length(dat2)]),2);
temp(1:length(dat1),1) = dat1;temp(1:length(dat2),2) = dat2;
colors = [0.75 0.75 0.75;[217 83 25]/255];
f;subplot(131),violinplot(temp,[],'ShowData',true,'ShowWhiskers',false,'ShowBox',false,'MarkerSize',5,'ViolinColor',colors);
set(gca,'tickdir','out','fontsize',16),box off; ylabel('Speed cm/s')

dat1 = Waves.baseline.ff.length;
dat2 = Waves.eOPN.ff.length;
temp = nan(max([length(dat1) length(dat2)]),2);
temp(1:length(dat1),1) = dat1;temp(1:length(dat2),2) = dat2;
colors = [0.75 0.75 0.75;[217 83 25]/255];
f;subplot(132),violinplot(temp,[],'ShowData',true,'ShowWhiskers',false,'ShowBox',false,'MarkerSize',5,'ViolinColor',colors);
set(gca,'tickdir','out','fontsize',16),box off;ylabel('Length \lambda')

dat1 = Waves.baseline.ff.duration;
dat2 = Waves.eOPN.ff.duration;
temp = nan(max([length(dat1) length(dat2)]),2);
temp(1:length(dat1),1) = dat1;temp(1:length(dat2),2) = dat2;
colors = [0.75 0.75 0.75;[217 83 25]/255];
f;subplot(133),violinplot(temp,[],'ShowData',true,'ShowWhiskers',false,'ShowBox',false,'MarkerSize',5,'ViolinColor',colors);
set(gca,'tickdir','out','fontsize',16),box off;ylabel('Duration ms')

f.Position = [681 559 760 420];
%%% Here we plot the wave directionality via histogram plot (looks nicer imo)
figure(2),histogram(Waves.baseline.ff.direction,-pi:pi/6:pi,'normalization','probability','edgecolor','none','facecolor',[0.5 0.5 0.5]),hold on
histogram(Waves.eOPN.ff.direction,-pi:pi/6:pi,'normalization','probability','edgecolor','none')
set(gca,'tickdir','out','fontsize',16),box off;ylabel('Probability');xlabel('Wave Direction')

%% FEEDBACK
f = figure(1);
dat1 = Waves.baseline.fb.speed;
dat2 = Waves.eOPN.fb.speed;
temp = nan(max([length(dat1) length(dat2)]),2);
temp(1:length(dat1),1) = dat1;temp(1:length(dat2),2) = dat2;
colors = [0.75 0.75 0.75;[217 83 25]/255];
f;subplot(131),violinplot(temp,[],'ShowData',true,'ShowWhiskers',false,'ShowBox',false,'MarkerSize',5,'ViolinColor',colors);
set(gca,'tickdir','out','fontsize',16),box off; ylabel('Speed cm/s')

dat1 = Waves.baseline.fb.length;
dat2 = Waves.eOPN.fb.length;
temp = nan(max([length(dat1) length(dat2)]),2);
temp(1:length(dat1),1) = dat1;temp(1:length(dat2),2) = dat2;
colors = [0.75 0.75 0.75;[217 83 25]/255];
f;subplot(132),violinplot(temp,[],'ShowData',true,'ShowWhiskers',false,'ShowBox',false,'MarkerSize',5,'ViolinColor',colors);
set(gca,'tickdir','out','fontsize',16),box off;ylabel('Length \lambda')

dat1 = Waves.baseline.fb.duration;
dat2 = Waves.eOPN.fb.duration;
temp = nan(max([length(dat1) length(dat2)]),2);
temp(1:length(dat1),1) = dat1;temp(1:length(dat2),2) = dat2;
colors = [0.75 0.75 0.75;[217 83 25]/255];
f;subplot(133),violinplot(temp,[],'ShowData',true,'ShowWhiskers',false,'ShowBox',false,'MarkerSize',5,'ViolinColor',colors);
set(gca,'tickdir','out','fontsize',16),box off;ylabel('Duration ms')

f.Position = [681 559 760 420];
%%% Here we plot the wave directionality via histogram plot (looks nicer imo)
figure(2),histogram(Waves.baseline.fb.direction,-pi:pi/6:pi,'normalization','probability','edgecolor','none','facecolor',[0.5 0.5 0.5]),hold on
histogram(Waves.eOPN.fb.direction,-pi:pi/6:pi,'normalization','probability','edgecolor','none')
set(gca,'tickdir','out','fontsize',16),box off;ylabel('Probability');xlabel('Wave Direction')
%% STATS table
stats.ff.speed = ranksum(Waves.baseline.ff.speed,Waves.eOPN.ff.speed);
stats.ff.length = ranksum(Waves.baseline.ff.length,Waves.eOPN.ff.length);
stats.ff.duration = ranksum(Waves.baseline.ff.duration,Waves.eOPN.ff.duration);
[stats.ff.direction, ~, ~] = circ_kuipertest(Waves.baseline.ff.direction, Waves.eOPN.ff.direction, 60, 0);


stats.fb.speed = ranksum(Waves.baseline.fb.speed,Waves.eOPN.fb.speed);
stats.fb.length = ranksum(Waves.baseline.fb.length,Waves.eOPN.fb.length);
stats.fb.duration = ranksum(Waves.baseline.fb.duration,Waves.eOPN.fb.duration);
[stats.fb.direction, ~, ~] = circ_kuipertest(Waves.baseline.fb.direction, Waves.eOPN.fb.direction, 60, 0);


%% Now look at LFP wave forms
if ~exist('LFP','var')
load('Z:\data\Austin\GridData\240508_eOPN_Inhib\LFPCombined.mat')
%load('Z:\data\Austin\GridData\240523_eOPNInhib\WavesCombined.mat')
end
addpath(genpath('C:\Users\khan332\Documents\GitHub\generalized-phase'));


[xgpbaseline,xobaseline] = lfpProcessing(LFP.baseline);
[xgpeOPN,xoeOPN] = lfpProcessing(LFP.eOPN);
%%
[PA.baseline.PA,PA.baseline.PA_angle] = calPhaseAlignment(xgpbaseline);

[PA.eOPN.PA,PA.eOPN.PA_angle] = calPhaseAlignment(xgpeOPN);

%%
d1 = xobaseline;
f = figure;
subplot(411),plot(squeeze(mean(mean(d1,3)))');xlim([100 1000]),axis off,hold on
subplot(4,1,2:4),imagesc(mean(d1,3)),colormap(redblue),caxis([-500 500]),xlim([100 1000])
hold on
for n = 1:32
    dat = squeeze(mean(d1(n,:,:),3));
plot(1:sz(3),n+(dat/max(dat)),'w'),hold on
end
f.Position = [681 159 360 720];
%%
d1 = PA.baseline.PA;
figure,imagesc(squeeze(d1));
figure,plot(squeeze(mean(d1,1)))
%% local functions
function output = remakeDat(dat)
% Make dat cells into labelled wave structures
% (start, source, speed, length, duration, direction, and evaluation
% point)
%ignore this
% output.start = dat{1};
% output.source = dat{2};
% output.speed = dat{3};
% output.length = dat{4};
% output.duration = dat{5};
% output.direction = dat{6};
% output.evaluation  = dat{7};
output.evaluation = dat{1};
output.source = dat{2};
output.start = dat{3};
output.duration = dat{4};
output.speed = dat{5};
output.direction = dat{6};
output.length  = dat{7};
end

function [xgp,xo] = lfpProcessing(dat)
t = permute(dat, [2 3 1]);
t = reshape(t,32,[]);
xo = bandpass_filter(t,5,40,1000); %x,f1,f2,Fs
sz = size(dat);
xo = reshape(xo,sz(2),1,sz(3)*sz(1));
[xgp,wt] = generalized_phase(xo, 1000, 0 );
xgp= reshape(xgp,sz(2),1,sz(3),sz(1));
xo= reshape(xo,sz(2),sz(3),sz(1));
temp = {};
for n = 1:105
    temp{n} = xgp(:,:,:,n);
end
xgp = temp;
end


