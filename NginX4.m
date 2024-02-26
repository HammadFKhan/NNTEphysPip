% clear; clc; 
% close all;
addpath(genpath('Main'));
% addpath(genpath('chronux'));
% addpath(genpath('Kilosort'));
addpath(genpath('npy-matlab'));
addpath(genpath('spikes-master'));
% IntanConcatenate legacy version
ds_filename = intanPreprocessing2(1); %IntanDs flag  %% double check file type
%% Combine intan data if needed
fpath = kilosortbinCombine();
%% Run Kilosort3 
% load only neccessary variables from memory mapped file
data = matfile(ds_filename);
fpath = data.fpath;
% Kilosort264FTestcode
% Kilosort264SharpTestcode
savepath = fullfile(fpath,['loadme','.mat']);
save(savepath,'ds_filename');
clearvars -except ds_filename
%% New load me
[fname,fpath] = uigetfile();
savepath = fullfile(fpath,['loadme','.mat']);
ds_filename = fullfile(fpath,fname);
save(savepath,'ds_filename');
%% Parameters for behaviour
data = matfile(ds_filename); % ds_filename comes from loadme.mat
% check if data directory matches where the file originated; if not we note
% the new directory path
parameters.experiment = 'cue'; % self - internally generated, cue - cue initiated
parameters.opto = 0; % 1 - opto ON , 0 - opto OFF
parameters.cool = 0; % No Cool 
parameters.windowBeforePull = 1.5; % in seconds
parameters.windowAfterPull = 1.5; % in seconds
parameters.windowBeforeCue = 1.5; % in seconds
parameters.windowAfterCue = 1.5; % in seconds
parameters.windowBeforeMI = 1.5; % in seconds 
parameters.windowAfterMI = 1.5; % in seconds 
parameters.Fs = 1000; % Eventual downsampled data
parameters.ts = 1/parameters.Fs;
parameters.IntanFs = data.targetedFs;
[Behaviour] = readLever(parameters,data.amplifierTime);
[IntanBehaviour] = readLeverIntan(parameters,data.amplifierTime,data.analogChannels(1,:),data.digitalChannels,Behaviour,1);
% Calculate ITI time for trials and reward/no reward sequence
temp1 = arrayfun(@(x) x.LFPtime(1), IntanBehaviour.cueHitTrace);
temp1 = vertcat(temp1,ones(1,IntanBehaviour.nCueHit)); %  write 1 for reward given
temp2 = arrayfun(@(x) x.LFPtime(1), IntanBehaviour.cueMissTrace);
temp2 = vertcat(temp2,zeros(1,IntanBehaviour.nCueMiss)); %  write 0 for no reward given
temp = [temp1,temp2];
[~,idx] = sort(temp(1,:)); %sort by occurance
IntanBehaviour.ITI = temp(:,idx);
%% For Cooling 
temperature = data.analogChannels(1,:);
temperature = (temperature-1.25)/0.005;
IntanBehaviour.temperature = resample(temperature,parameters.Fs,data.targetedFs);
clear temperature
for n = 1:IntanBehaviour.nCueHit
    hitTemp(n,1) = IntanBehaviour.temperature(IntanBehaviour.cueHitTrace(n).LFPIndex(1));
end
for n = 1:IntanBehaviour.nCueMiss
    missTemp(n,1) = IntanBehaviour.temperature(IntanBehaviour.cueMissTrace(n).LFPIndex(1));
end
for n = 1:length(IntanBehaviour.missTrace)
    FATemp(n,1) = IntanBehaviour.temperature(IntanBehaviour.missTrace(n).LFPIndex(1));
end
%% Plot behaviour
figure
for i=1:IntanBehaviour.nCueHit
    plot(0:3000,smoothdata(IntanBehaviour.cueHitTrace(i).trace),'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
    try
        hitTrace(i,:) = smoothdata(IntanBehaviour.cueHitTrace(i).rawtrace);
    catch
        continue
    end
end
for n = 1:IntanBehaviour.nCueHit
    IntanBehaviour.AvgHitTrace(n,:) = IntanBehaviour.cueHitTrace(n).trace;
end
IntanBehaviour.AvgHitTrace = mean(IntanBehaviour.AvgHitTrace,1);
for n = 1:IntanBehaviour.nCueMiss
    IntanBehaviour.AvgMissTrace(n,:) = IntanBehaviour.cueMissTrace(n).trace;
end
IntanBehaviour.AvgMissTrace = mean(IntanBehaviour.AvgMissTrace,1);
IntanBehaviour.AvgHitTrace = mean(IntanBehaviour.AvgHitTrace,1);

%% LFP probe setup for 64F and analysis
% Since there are two probes we want to seperate everything into linear
% maps for CSD and depthwise LFP analysis and then we do filtering
data = matfile(ds_filename);
% load UCLA_chanMap_64F2
load UCLA_chanmap_fixed.mat
if ~exist('lfp','var'),lfp = data.amplifierData;end
%TODO check if the field orientation during insertion is reversed (ie. probe 1 is lateral to probe 2)
probe1 = lfp(s.sorted_probe_wiring(:,5)==1,:);
probe2 = lfp(s.sorted_probe_wiring(:,5)==2,:);
chanProbe1 = s.sorted_probe_wiring(s.sorted_probe_wiring(:,5)==1,:); %needed for linear channel mapping later
chanProbe2 = s.sorted_probe_wiring(s.sorted_probe_wiring(:,5)==2,:); 
clear lfp
% LFP filter
set(0,'DefaultFigureWindowStyle','normal')
LFP.probe1= fastpreprocess_filtering(probe1,data.targetedFs);
% LFP.probe1 = bestLFP(LFP.probe1);
% LFP.probe1 = bandFilter(LFP.probe1,'depth'); % Extract LFPs based on 'depth' or 'single'
if ~isempty(probe2)
    LFP.probe2 = fastpreprocess_filtering(probe2,data.targetedFs);
    LFP.probe2 = bestLFP(LFP.probe2);
    LFP.probe2 = bandFilter(LFP.probe2,'depth'); % Extract LFPs based on 'depth' or 'single'
end
%%% Calculate generalize phase for electrodes
addpath(genpath('C:\Users\khan332\Documents\GitHub\generalized-phase'));
xo = bandpass_filter(LFP.probe1.LFP,5,40,1000); %x,f1,f2,Fs
sz = size(xo);
xo = reshape(xo,sz(1),1,sz(2));
linearProbe = find(s.sorted_probe_wiring(:,2)==20);
[xgp,wt] = generalized_phase(xo, 1000, 0 );
LFP.probe1.xgp = xgp;
LFP.probe1.xo = xo;
%%% Save data
[fpath,name,exts] = fileparts(ds_filename);
sessionName = [fpath,'/','LFP.mat'];
% save(sessionName,"IntanBehaviour","fpath","parameters","-v7.3");
save(sessionName,"IntanBehaviour","parameters","LFP","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
%%
for i = 1:IntanBehaviour.nCueHit
    disp(['Analyzing trial: ' num2str(i)])
    hitWin = floor(IntanBehaviour.cueHitTrace(i).LFPtime*1000); %multiply by Fs;
    hitxgp{i} = LFP.probe1.xgp(:,hitWin);
    hitLFP{i} = LFP.probe1.xo(:,hitWin);
end
for i = 1:IntanBehaviour.nCueMiss
    disp(['Analyzing trial: ' num2str(i)])
    missWin = floor(IntanBehaviour.cueMissTrace(i).LFPtime*1000);
    missxgp{i} = LFP.probe1.xgp(:,missWin);
    missLFP{i} = LFP.probe1.xo(:,missWin);
end

%%% Motion Intiated
for i = 1:length(IntanBehaviour.MIHitTrace)
    disp(['Analyzing trial: ' num2str(i)])
    MIhitWin = floor(IntanBehaviour.MIHitTrace(i).LFPtime*1000); %multiply by Fs;
    MIhitxgp{i} = LFP.probe1.xgp(:,MIhitWin);
    MIhitLFP{i} = LFP.probe1.xo(:,MIhitWin);
end

for i = 1:length(IntanBehaviour.MIFATrace)
    disp(['Analyzing trial: ' num2str(i)])
    MIFAWin = floor(IntanBehaviour.MIFATrace(i).LFPtime*1000);
    MIFAxgp{i} = LFP.probe1.xgp(:,MIFAWin);
    MIFALFP{i} = LFP.probe1.xo(:,MIFAWin);
end
hitxgp = cellfun(@(x) reshape(x,64,1,[]),hitxgp,'UniformOutput',false);
missxgp = cellfun(@(x) reshape(x,64,1,[]),missxgp,'UniformOutput',false);
MIhitxgp = cellfun(@(x) reshape(x,64,1,[]),MIhitxgp,'UniformOutput',false);
MIFAxgp = cellfun(@(x) reshape(x,64,1,[]),MIFAxgp,'UniformOutput',false);
%%
load myMap
[hitPA,hitPA_angle] = calPhaseAlignment(hitxgp);
[missPA,missPA_angle] = calPhaseAlignment(missxgp);
[MIhitPA,MIhitPA_angle] = calPhaseAlignment(MIhitxgp);
[MIFAPA,MIFAPA_angle] = calPhaseAlignment(MIFAxgp);
% figure,plot((squeeze(PA)'),'LineWidth',1);
% cmap = (gray(22));
% set(gca(),'ColorOrder',cmap)
% cmap = colormap_redblackblue();
%% PA plotting
linearProbe = find(s.sorted_probe_wiring(:,2)==20);
figure,
subplot(311),imagesc(-1500:1500,1:22,squeeze((hitPA(linearProbe,:,1:3001)))),xlim([-500 1500]),colormap(myMap),colorbar,title('Hit')
subplot(312),imagesc(-1500:1500,1:22,squeeze((missPA(linearProbe,:,1:3001)))),xlim([-500 1500]),colormap(myMap),colorbar,title('Miss')
% subplot(313),imagesc(-1500:1500,1:22,squeeze((MIhitPA(linearProbe,:,1:3001)))),xlim([-500 1500]),colormap(myMap),colorbar,title('MI Hit')
subplot(313),imagesc(-1500:1500,1:22,squeeze((MIFAPA(linearProbe,:,1:3001)))),xlim([-500 1500]),colormap(myMap),colorbar,title('MI FA')
%% Alternative way of ploting PA w/angle
figure,
shuf = hitPA_angle(:);
t = squeeze(hitPA_angle);
subplot(331),histogram(t,-pi:pi/8:pi,'normalization','probability'),set(gca,'TickDir','out'),set(gca,'fontsize',12),box off
t = squeeze(missPA_angle);
subplot(334),histogram(t,-pi:pi/8:pi,'normalization','probability'),set(gca,'TickDir','out'),set(gca,'fontsize',12),box off
[~,r] = max(squeeze(MIFAPA),[],2);
t = squeeze(MIFAPA_angle);
subplot(337),histogram(t,-pi:pi/8:pi,'normalization','probability'),set(gca,'TickDir','out'),set(gca,'fontsize',12),box off
sz = size(squeeze(hitPA));
dat = mean(squeeze(hitPA(linearProbe,:,1:3001)));
subplot(3,3,[2 3]),plot(1:sz(2),dat);set(gca,'TickDir','out'),set(gca,'fontsize',12),box off
dat = mean(squeeze(missPA(linearProbe,:,1:3001)));
subplot(3,3,[5 6]),plot(1:sz(2),dat);set(gca,'TickDir','out'),set(gca,'fontsize',12),box off
dat = mean(squeeze(MIFAPA(linearProbe,:,1:3001)));
subplot(3,3,[8 9]),plot(1:sz(2),dat);set(gca,'TickDir','out'),set(gca,'fontsize',12),box off
%%
figure,
t = squeeze(missPA_angle);
dat = mean(smoothdata(t(linearProbe,:),'gaussian',5),2);
err = std(t(linearProbe,:),[],2)/sqrt(size(t,2));
subplot(3,2,[1 3 5]),errorbar(dat,1:length(linearProbe),err,'.-','horizontal'),axis tight,set(gca,'TickDir','out'),set(gca,'fontsize',12),box off
set(gca, 'YDir','reverse')
subplot(322),histogram(t(linearProbe(4),:),-pi:pi/8:pi,'normalization','probability'),set(gca,'TickDir','out'),set(gca,'fontsize',12),box off,title('Superficial')
subplot(324),histogram(t(linearProbe(10),:),-pi:pi/8:pi,'normalization','probability'),set(gca,'TickDir','out'),set(gca,'fontsize',12),box off,title('Middle')
subplot(326),histogram(t(linearProbe(20),:),-pi:pi/8:pi,'normalization','probability'),set(gca,'TickDir','out'),set(gca,'fontsize',12),box off,title('Deep')
%% Spikes analysis
[fpath,name,exts] = fileparts(ds_filename);
data = matfile(ds_filename);
path = [fpath,'/kilosort3/'];
mergename = 'merged';
Kilosort3AutoMergeTester
path = [fpath,'/kilosort3/' mergename];
%%% Spike preprocessing (includes merging (optional) and channel info
%%% return)
% Read in kilosort data for matlab analysis
SpikeClusters = readNPY(fullfile(path, 'spike_clusters.npy'));
SpikeSamples = readNPY(fullfile(path, 'spike_times.npy'));
SpikeChannel = readNPY(fullfile(path,'channel_positions.npy'));
Spikes.SpikeClusters = SpikeClusters; 
Spikes.SpikeSamples = SpikeSamples;
Spikes = clusterSort(Spikes);
Spikes = ISI(Spikes,0.01,data.Fs,0); %Spikes, Interval, Fs
% Calculate Depth profile
% load chanMap64F2
load UCLA_chanMap
[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms, max_site] =...
    spikeTemplatePosition(data.fpath,ycoords,[]);
for i = 1:length(tempAmps)
    Spikes.Clusters(i).spikeDepth = templateDepths(i);
    Spikes.Clusters(i).channelDepth = max_site(i);
    Spikes.Clusters(i).spikeAmplitude = tempAmps(i);
    Spikes.Clusters(i).waveforms = waveforms(i,:);
    Spikes.Clusters(i).spikeDuration = templateDuration(i)/data.Fs*1000;
end
%%% delete empty spikes
temp = arrayfun(@(x) isempty(x.cluster), Spikes.Clusters);
Spikes.Clusters(temp) = []; 
%% Calculate trial PSTH for lever
Spikes = leverPSTH(Spikes,IntanBehaviour);
%%% save spike output data to load into gui
savepath = fullfile(path,['spks4sorting','.mat']);
save(savepath,'Spikes')
ManualSpikeCurateGUI
%% Basic spike analysis
% z-score spike rates
if exist('goodSpkComponents','var')
    Spikes.goodSpkComponents = unique(goodSpkComponents);
end
[Spikes] = sortSpkLever(Spikes,IntanBehaviour1);
[fpath,name,exts] = fileparts(ds_filename);
sessionName = [fpath,'/','Spikes.mat'];
% save(sessionName,"IntanBehaviour","fpath","parameters","-v7.3");
save(sessionName,"Spikes","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
%% Make it layer specific TODO:Gamma GED/Spike coherence
Spikes = layerspikeAnalysis(Spikes,IntanBehaviour,LFP);
%%
figure,hold on
for n = 20
subplot(2,1,[1]),Show_Spikes(Spikes.PSTH.hit.spks{n}),axis off
subplot(2,1,[2]),bar(-1500:1500,smoothdata(Spikes.PSTH.hit.spkRates(n,:)),'FaceColor',[28/255 117/255 188/255],'EdgeColor','none')
axis tight, box off, set(gca,'TickDir','out')
set(gca,'fontsize',16)
xline(0)
xline(mean(IntanBehaviour.reactionTime)*1000)
end
%%
figure,subplot(121),imagesc(-IntanBehaviour.parameters.windowBeforeCue*1000:IntanBehaviour.parameters.windowAfterCue*1000,...
    1:size(Spikes.PSTH.hit.normSpk(Spikes.PSTH.hit.spkIdx),1),Spikes.PSTH.hit.normSpk(Spikes.PSTH.hit.spkIdx,:)),hold on
colormap(flip(gray))
colorbar
set(gca,'fontsize',16)
caxis([0.0 2.56])
subplot(122),imagesc(-IntanBehaviour.parameters.windowBeforeCue*1000:IntanBehaviour.parameters.windowAfterCue*1000,...
    1:size(Spikes.PSTH.miss.normSpk(Spikes.PSTH.miss.spkIdx),1),Spikes.PSTH.miss.normSpk(Spikes.PSTH.hit.spkIdx,:)),hold on
colormap(flip(gray))
colorbar
set(gca,'fontsize',16)
caxis([0.0 2.56])
%% Neural Trajectory Segementation using GPFA
% Note that we concatenate trial conditions as to apply the same models for
% statistical comparison (ie. hit vs miss, hit vs FA, opto vs no opto)
Spikes = makeSpikeGPFA(Spikes);
Spikes.GPFA.HitMiss.dat = [Spikes.GPFA.hit.dat,Spikes.GPFA.miss.dat];
for n = IntanBehaviour.nCueHit+1:length(Spikes.GPFA.HitMiss.dat) %fix trials
    Spikes.GPFA.HitMiss.dat(n).trialId = n;
end
Spikes.GPFA.MIHitFA.dat = [Spikes.GPFA.MIHit.dat,Spikes.GPFA.MIFA.dat];
for n = length(IntanBehaviour.MIHitTrace)+1:length(Spikes.GPFA.MIHitFA.dat) %fix trials
    Spikes.GPFA.MIHitFA.dat(n).trialId = n;
end
%%
addpath(genpath('C:\Users\khan332\Documents\GitHub\NeuralTraj'));
addpath(genpath('mat_results'));
if exist('mat_results','dir'),rmdir('mat_results','s'),end
[Spikes.GPFA.resultHit,Spikes.GPFA.seqTrainHit] = gpfaAnalysis(Spikes.GPFA.hit.dat,1); %Run index
[Spikes.GPFA.resultMiss,Spikes.GPFA.seqTrainMiss] = gpfaAnalysis(Spikes.GPFA.miss.dat,2); %Run index
[Spikes.GPFA.resultMIHit,Spikes.GPFA.seqTrainMIHit] = gpfaAnalysis(Spikes.GPFA.MIHit.dat,3); %Run index
[Spikes.GPFA.resultMIFA,Spikes.GPFA.seqTrainMIFA] = gpfaAnalysis(Spikes.GPFA.MIFA.dat,4); %Run index
[Spikes.GPFA.resultHitMiss,Spikes.GPFA.seqTrainHitMiss] = gpfaAnalysis(Spikes.GPFA.HitMiss.dat,5); %Run index
[Spikes.GPFA.resultMIHitFA,Spikes.GPFA.seqTrainMIHitFA] = gpfaAnalysis(Spikes.GPFA.MIHitFA.dat,6); %Run index
%% Neural Trajectory Analysis
neuralTrajAnalysis(Spikes,[],IntanBehaviour)
%%
figure,
for n = 1:111
    subplot(131),plot(-74*20:20:20*75,Spikes.GPFA.seqTrainHit(n).xorth(1,:),'Color',[0 0 0 0.4]),hold on,box off,set(gca,'fontsize',18),xlim([-1500 1500]),ylim([-1 2])
    subplot(132),plot(-74*20:20:20*75,Spikes.GPFA.seqTrainHit(n).xorth(2,:),'Color',[0 0 0 0.4]),hold on,box off,set(gca,'fontsize',18),xlim([-1500 1500]),ylim([-1 2])
    subplot(133),plot(-74*20:20:20*75,Spikes.GPFA.seqTrainHit(n).xorth(3,:),'Color',[0 0 0 0.4]),hold on,box off,set(gca,'fontsize',18),xlim([-1500 1500]),ylim([-1 3])
end





%% ----------------------- Some local function to make things easier ----------------------- %%
% Hit trials
function output = leverLFPAnalysis(linearProbe,Behaviour) % LFP of linear channel and Behaviour struct
CSDoutputhit = [];waveletHit = [];waveletMiss = [];powerCWThit = [];CSDoutputmiss = []; hitLFP = [];missLFP = [];
params.tapers = [5 9];
params.Fs = 1000;
params.fpass = [4 80];
params.err = [2 0.05];
timpnts = size(Behaviour.AvgHitTrace,2); %timepoints
for i = 1:Behaviour.nCueHit
    disp(['Analyzing trial: ' num2str(i)])
    timestamphit(i,:) = [Behaviour.cueHitTrace(i).LFPtime(1),Behaviour.cueHitTrace(i).LFPtime(end)]; %taken in seconds
    hitWin = floor(Behaviour.cueHitTrace(i).LFPtime*1000); %multiply by Fs;
    hitLFP(:,:,i) = linearProbe(:,hitWin);
    for n = 1:size(hitLFP,1)
        [temp(:,:,n), fwt] = calCWTSpectogram(hitLFP(n,:,i),1:timpnts,1000,10,[10 40],0,1);
    end
    powerCWThit(:,:,i) = mean(temp,3);
    [CSDoutputhit(:,:,i)]  = CSD(hitLFP(:,:,i)'/1E6,1000,50E-6);
end
for i = 1:Behaviour.nCueMiss
    disp(['Analyzing trial: ' num2str(i)])
    timestampmiss(i,:) = [Behaviour.cueMissTrace(i).LFPtime(1),Behaviour.cueMissTrace(i).LFPtime(end)]; %taken in seconds
    missWin = floor(Behaviour.cueMissTrace(i).LFPtime*1000);
    missLFP(:,:,i) = linearProbe(:,missWin);
    for n = 1:size(missLFP,1)
        [temp(:,:,n), fwt] = calCWTSpectogram(missLFP(n,:,i),1:timpnts,1000,10,[10 40],0,1);
    end
    powerCWTmiss(:,:,i) = mean(temp,3);
    [CSDoutputmiss(:,:,i)]  = CSD(missLFP(:,:,i)'/1E6,1000,20E-6);
end
%%% Motion Initiated Analysis for FA
for i = 1:length(Behaviour.MIHitTrace)
    disp(['Analyzing trial: ' num2str(i)])
    timestampMIHit(i,:) = [Behaviour.MIHitTrace(i).LFPtime(1),Behaviour.MIHitTrace(i).LFPtime(end)]; %taken in seconds
    MIHitWin = floor(Behaviour.MIHitTrace(i).LFPtime*1000);
    MIHitLFP(:,:,i) = linearProbe(:,MIHitWin);
    for n = 1:size(MIHitLFP,1)
        [temp(:,:,n), fwt] = calCWTSpectogram(MIHitLFP(n,:,i),1:timpnts,1000,10,[10 40],0,1);
    end
    powerCWTMIHit(:,:,i) = mean(temp,3);
    [CSDoutputMIHit(:,:,i)]  = CSD(MIHitLFP(:,:,i)'/1E6,1000,20E-6);

end
for i = 1:length(Behaviour.MIFATrace)
    disp(['Analyzing trial: ' num2str(i)])
    timestampMIFA(i,:) = [Behaviour.MIFATrace(i).LFPtime(1),Behaviour.MIFATrace(i).LFPtime(end)]; %taken in seconds
    MIFAWin = floor(Behaviour.MIFATrace(i).LFPtime*1000);
    MIFALFP(:,:,i) = linearProbe(:,MIFAWin);
    for n = 1:size(missLFP,1)
        [temp(:,:,n), fwt] = calCWTSpectogram(MIFALFP(n,:,i),1:timpnts,1000,10,[10 40],0,1);
    end
    powerCWTMIFA(:,:,i) = mean(temp,3);
    [CSDoutputMIFA(:,:,i)]  = CSD(MIFALFP(:,:,i)'/1E6,1000,20E-6);
end
for i = 1:size(linearProbe,1)
    disp(['Analyzing electrode: ' num2str(i)])
    tfmiss(:,:,i) = itpc(linearProbe(i,:),timestampmiss,1000,0);
    [tfhit(:,:,i),frex,pnts] = itpc(linearProbe(i,:),timestamphit,1000,0);
    [tfMIhit(:,:,i),frex,pnts] = itpc(linearProbe(i,:),timestampMIHit,1000,0);
    [tfMIFA(:,:,i),frex,pnts] = itpc(linearProbe(i,:),timestampMIFA,1000,0);
end
% figure(), clf
% contourf(1:pnts,frex,mean(tfhit,3),10,'linecolor','none')
% set(gca,'clim',[0 .2],'ydir','normal')
% title('ITPC')
% colormap(jet),colorbar
% figure,plot(smoothdata(squeeze(mean(tfhit,1)),'gaussian',50),'k')
% hold on,plot(smoothdata(squeeze(mean(tfmiss,1)),'gaussian',50),'r')
[Shit,fhit,Serrhit]=mtspectrumc(squeeze(mean(hitLFP,1)),params);
[Smiss,~,Serrmiss] = mtspectrumc(squeeze(mean(missLFP,1)),params);

% output data
output.hitLFP = hitLFP;
output.powerCWThit = powerCWThit;
output.tfhit = tfhit;
output.fwt = fwt;
output.frex = frex;
output.pnts = pnts;
output.Shit = Shit;

output.Sf = fhit;
output.Serrhit = Serrhit;
output.missLFP = missLFP;
output.powerCWTmiss = powerCWTmiss;
output.tfmiss = tfmiss;
output.Smiss = Smiss;
output.Serrmiss = Serrmiss;

output.MIHitLFP = MIHitLFP;
output.powerCWTMIHit = powerCWTMIHit;
output.tfMIhit = tfMIhit;

output.MIFALFP = MIFALFP;
output.powerCWTMIFA = powerCWTMIFA;
output.tfMIFA = tfMIFA;

end

function [ampHit,ampMiss] = calcBandPower(hitBand,missBand)
narrowBand = [hitBand,missBand];
n_iter = 1000;
nhits = size(hitBand,2);
ntotal = size(narrowBand,2);
nullDistHit = zeros(size(narrowBand,1),n_iter);
nullDistMiss = zeros(size(narrowBand,1),n_iter);
for n = 1:ntotal %trial
    temp  = narrowBand(:,n);
    amp(:,n) = abs(hilbert(temp)).^2;
end
ampHit = amp(:,1:nhits);
ampMiss = amp(:,nhits+1:end);
fprintf('Shuffling with %d interations...',n_iter);
for n = 1:n_iter
    randIdx = randperm(ntotal);
    hitRand = amp(:,randIdx(1:nhits));
    missRand = amp(:,randIdx(nhits+1:end));
    nullDistHit(:,n) = mean(hitRand,2);
    nullDistMiss(:,n) = mean(missRand,2);
end
fprintf('done\n');
% now use permutation testing to get Z-value
% the value we use is the normalized distance away from the mean of
% boot-strapped values
ampHit = (ampHit-mean(nullDistHit,2))./std(nullDistHit,[],2);
ampMiss = (ampMiss-mean(nullDistMiss,2))./std(nullDistMiss,[],2);
ampHit = ampHit';
ampMiss = ampMiss';
end

function [prestim,poststim,prestimIdx,poststimIdx] = stimulusPower(data,Behaviour)
for n = 1:size(data,1) %trials
    % Prestimulus
    win = 1:Behaviour.parameters.windowBeforeCue*1000;
    for nn = 1:size(data,3) %electrodes
        prestim(n,nn) = squeeze(mean(data(n,win,nn),2));
        [~,prestimIdx(n,nn)] = max(data(n,win,nn));
    end
    % post stimulus
    win = 1+Behaviour.parameters.windowBeforeCue*1000:Behaviour.parameters.windowAfterCue*1000;
    for nn = 1:size(data,3)
        poststim(n,nn) = squeeze(mean(data(n,win,nn),2));
        [~,poststimIdx(n,nn)] = max(data(n,win,nn));
    end
end
end
function output = GEDITPC(GED)
broadBandIdx = 26:54; %Broadband idx base on spectrogram.frex (35:40),30:35 low Beta, 55:60 gamma
temp = GED.comp1;
%Search freq space and note maximal response for beta
betatfhit = temp.tfhit(broadBandIdx,:);
[~,idx] = max(betatfhit(:,:),[],1);
[n,bin] = hist(idx,unique(idx));
[~,idx] = sort(-n);
tFreq = broadBandIdx(bin(idx(1:10)));
output.comp1.hit = smoothdata(squeeze(temp.tfhit(tFreq,:)),2);
output.comp1.miss = smoothdata(squeeze(temp.tfmiss(tFreq,:)),2);
output.comp1.MIhit = smoothdata(squeeze(temp.tfMIhit(tFreq,:)),2);
output.comp1.MIFA = smoothdata(squeeze(temp.tfMIFA(tFreq,:)),2);

temp = GED.comp2;
betatfhit = temp.tfhit(broadBandIdx,:);
[~,idx] = max(betatfhit(:,:),[],1);
[n,bin] = hist(idx,unique(idx));
[~,idx] = sort(-n);
tFreq = broadBandIdx(bin(idx(1:10)));
output.comp2.hit = smoothdata(squeeze(temp.tfhit(tFreq,:)),2);
output.comp2.miss = smoothdata(squeeze(temp.tfmiss(tFreq,:)),2);
output.comp2.MIhit = smoothdata(squeeze(temp.tfMIhit(tFreq,:)),2);
output.comp2.MIFA = smoothdata(squeeze(temp.tfMIFA(tFreq,:)),2);

temp = GED.comp3;
betatfhit = temp.tfhit(broadBandIdx,:);
[~,idx] = max(betatfhit(:,:),[],1);
[n,bin] = hist(idx,unique(idx));
[~,idx] = sort(-n);
tFreq = broadBandIdx(bin(idx(1:10)));
output.comp3.hit = smoothdata(squeeze(temp.tfhit(tFreq,:)),2);
output.comp3.miss = smoothdata(squeeze(temp.tfmiss(tFreq,:)),2);
output.comp3.MIhit = smoothdata(squeeze(temp.tfMIhit(tFreq,:)),2);
output.comp3.MIFA = smoothdata(squeeze(temp.tfMIFA(tFreq,:)),2);
end

function [trials,win] = makeSpikeWin(Spikes,spikeId,Fs)
%% Begin analysis across clusters
winLen = 4;
win = double(0:winLen:max(Spikes.SpikeSamples)/Fs);
trials = cell(1,length(win));
% win(1:400) = [];
for i = 1:length(win)-1
    spiketm = Spikes.Clusters(spikeId).spikeTime-win(i);
    idx = find(spiketm>=0 & spiketm < winLen);
    trials{i} = zeros(1,length(0:0.001:winLen));
    spiked = discretize(spiketm(idx),0:0.001:winLen);
    if ~isnan(spiked)
        trials{i}(spiked) = 1;
    end
end
end





function plotspikeCoherence(spikea,spikeb,spikeCoherence)

interval = 0:2000;
spiker = spikea;

figure;
%****************
disp('Plotting rastergrams (slow) ...');
subplot('position',[0.1 0.4 0.4 0.55]);
make_nice_spike_raster(spiker);
V = axis;
axis([0 size(spikea{1},2) V(3) V(4)]);
grid on;
ylabel('Trial Number');
title(fprintf('Unit rasters'));
%*************
subplot('position',[0.1 0.1 0.4 0.3]);
smooth_window = 25;  % give sigma of 12.5ms
make_nice_mean_raster(spiker,smooth_window);
V = axis;
axis([0 size(spikea{1},2) V(3) V(4)]);
plot([interval(1),interval(1)],[V(3),V(4)],'k-'); hold on;
plot([interval(end),interval(end)],[V(3),V(4)],'k-'); hold on;
ylabel('Firing Rate');
xlabel('Time (ms)');

spiker = spikeb;
%****************
disp('Plotting rasters2 per trial (slow) ...');
subplot('position',[0.57 0.4 0.4 0.55]);
make_nice_spike_raster(spiker);
V = axis;
axis([0 size(spikeb{1},2) V(3) V(4)]);
grid on;
ylabel('Trial Number');
title(fprintf('Unit rasters'));
%*************
subplot('position',[0.57 0.1 0.4 0.3]);
smooth_window = 25;  % give sigma of 12.5ms
make_nice_mean_raster(spiker,smooth_window);
V = axis;
axis([0 size(spikeb{1},2) V(3) V(4)]);
plot([interval(1),interval(1)],[V(3),V(4)],'k-'); hold on;
plot([interval(end),interval(end)],[V(3),V(4)],'k-'); hold on;
ylabel('Firing Rate');
xlabel('Time (ms)');

figure;
colo = 'rbgy';
CNUM=1;
subplot(2,2,1);
for ii = 1:CNUM
    H = semilogx(spikeCoherence.freq{ii},spikeCoherence.coho{ii},[colo(ii),'-']); hold on;
    set(H,'Linewidth',2);
    H = semilogx(spikeCoherence.freq{ii},(spikeCoherence.coho{ii}+spikeCoherence.scoho{ii}),[colo(ii),':']); hold on;
    H = semilogx(spikeCoherence.freq{ii},(spikeCoherence.coho{ii}-spikeCoherence.scoho{ii}),[colo(ii),':']);
    H = semilogx(spikeCoherence.freq{ii},spikeCoherence.rcoho{ii},[colo(ii),'--']); hold on;
    set(H,'Linewidth',1);
end
ylabel('Coherence Magnitude');
xlabel('Frequency (Hz)');
title(sprintf('Spike A with Spike B'));

subplot(2,2,2);
for ii = 1:CNUM
    H = plot(spikeCoherence.freq{ii},spikeCoherence.phaso{ii},[colo(ii),'-']); hold on;
    set(H,'Linewidth',2);
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.phaso{ii}+spikeCoherence.sphaso{ii}),[colo(ii),':']); hold on;
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.phaso{ii}-spikeCoherence.sphaso{ii}),[colo(ii),':']);
    %H = plot(freq{ii},rphaso{ii},[colo(ii),'--']); hold on;
    set(H,'Linewidth',1);
end
ylabel('Coherence Angle');
xlabel('Freq');
title(sprintf('Spike A with Spike B'));

subplot(2,2,3);
for ii = 1:CNUM
    H = plot(spikeCoherence.freq{ii},spikeCoherence.spika_pow{ii},[colo(ii),'-']); hold on;
    set(H,'Linewidth',2);
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.spika_pow{ii}+spikeCoherence.sspika_pow{ii}),[colo(ii),':']); hold on;
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.spika_pow{ii}-spikeCoherence.sspika_pow{ii}),[colo(ii),':']);
end
ylabel('Spike Unit A Pow');
xlabel('Frequency (Hz)');
title(sprintf('Spike A with Spike B'));

subplot(2,2,4);
for ii = 1:CNUM
    H = plot(spikeCoherence.freq{ii},spikeCoherence.spike_pow{ii},[colo(ii),'-']); hold on;
    set(H,'Linewidth',2);
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.spike_pow{ii}+spikeCoherence.sspike_pow{ii}),[colo(ii),':']); hold on;
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.spike_pow{ii}-spikeCoherence.sspike_pow{ii}),[colo(ii),':']);
end
ylabel('Spike Unit B Pow');
xlabel('Frequency (Hz)');
title(sprintf('Spike B with Spike A'));

end
%****************** function to make a nice raster plot ****************
function smorate = make_nice_mean_raster(spmat,smooth_window)
%*********** spmat1 and spmat2 are spike matrices of two conditions you wish to compare
%*********** smooth_window ... gaussian smoothing in millisecs
numconds = size(spmat,2);
if (numconds==2)
    colo = [[1,0,0];[0,0,1]];
else
    colo = [[1,0,0];[0,0,1];[0,1,0];[1,1,0]];
end

for k = 1:numconds
    spud = spmat{1,k};
    numtrials(k) = size(spud,1);
    smorate = gauss_smooth(sum( spud(1:numtrials(k),:))/....
        numtrials(k),smooth_window)*1000;
    [mSmo,idx] = max(smorate);
    if idx<1400 && idx>200
    smorate = gaussmf(1:length(smorate),[mSmo idx])'*mSmo;
    end
%     smorate = (smorate-min(smorate))/(max(smorate)-min(smorate));
%     smorate = zscore(smorate);
    H = plot(smorate,'k-'); hold on;
    set(H,'Color',colo(k,:));
end

end


%******************* make a rastergram of the actual spikes on each trial
function make_nice_spike_raster(spmat)

colo = [[1,1,1];[0,0,0]];
colormap(colo);

totspike = [];
for cubo = 1:size(spmat,2)
    totspike = [totspike; (cubo*spmat{1,cubo}) ];
    totspike = [totspike; (cubo*spmat{1,cubo}) ];
    
end
totspike = totspike + 1;
imagesc(totspike);
V = axis;
axis([V(1) V(2) 0 size(totspike,1)]);

end

%**************************************************************
function output = gauss_smooth(input, window)
% Smoothing function:
% output = smooth(input, window)
% "Window" is the total kernel width.
% Input array must be one-dimensional.

input_dims = ndims(input);
input_size = size(input);
if input_dims > 2 | min(input_size) > 1,
    disp('Input array is too large.');
    return
end

if input_size(2) > input_size(1),
    input = input';
    toggle_dims = 1;
else
    toggle_dims = 0;
end

if window/2 ~= round(window/2),
    window = window + 1;
end
halfwin = window/2;

input_length = length(input);
%********* gauss window +/- 1 sigma
x = -halfwin:1:halfwin;
kernel = exp(-x.^2/(window/2)^2);
kernel = kernel/sum(kernel);

padded(halfwin+1:input_length+halfwin) = input;
padded(1:halfwin) = ones(halfwin, 1)*input(1);
padded(length(padded)+1:length(padded)+halfwin) = ones(halfwin, 1)*input(input_length);

output = conv(padded, kernel);
output = output(window:input_length+window-1);

if toggle_dims == 1,
    output = output';
end
end
function plotBetaGammaPower(Behaviour,LFP)
figure,hold on
plotWin = -Behaviour.windowBeforeCue*1000:Behaviour.windowAfterCue*1000;
nChan = size(LFP.LFP,1);
for electrode = 1:nChan
    subplot(3,1,[1 2]),plot(plotWin,smoothdata((squeeze(mean(LFP.narrowband.hitLFP.betaPower(:,:,electrode),1))),'gaussian',10),'k','LineWidth',1),hold on
    plot(plotWin,smoothdata((squeeze(mean(LFP.narrowband.hitLFP.gammaPower(:,:,electrode),1))),'gaussian',10),'Color',[188/255 190/255 192/255],'LineWidth',1)
end
subplot(3,1,[1 2]),yline(2.58,'k--') % 99% cutoff
subplot(3,1,[1 2]),yline(-2.58,'k--') % 99% cutoff
subplot(3,1,[1 2]),yline(1.96,'k--') % 99% cutoff
subplot(3,1,[1 2]),yline(-1.96,'k--') % 99% cutoff
legend('Beta','Gamma')
box off,set(gca,'TickDir','out')
set(gca,'xtick',[])
set(gca,'FontSize',16)
ylabel('Power (Z)')
title('Hit Trials')
subplot(3,1,3),semilogy(plotWin,smoothdata(mean(LFP.narrowband.hitLFP.betaGammaPvalue,1),'gaussian',1),'k','LineWidth',1),hold on
box off,set(gca,'TickDir','out')
set(gca,'FontSize',16)
xlabel('Time from cue onset (ms)')
ylabel('P-value')

figure,hold on
for electrode = 1:nChan
subplot(3,1,[1 2]),plot(plotWin,smoothdata(mean(LFP.narrowband.missLFP.betaPower(:,:,electrode),1),'gaussian',10),'k','LineWidth',1),hold on
plot(plotWin,smoothdata(mean(LFP.narrowband.missLFP.gammaPower(:,:,electrode),1),'gaussian',10),'Color',[188/255 190/255 192/255],'LineWidth',1)
end
subplot(3,1,[1 2]),yline(2.58,'k--') % 99% cutoff
subplot(3,1,[1 2]),yline(-2.58,'k--') % 99% cutoff
subplot(3,1,[1 2]),yline(1.96,'k--') % 99% cutoff
subplot(3,1,[1 2]),yline(-1.96,'k--') % 99% cutoff
legend('Beta','Gamma')
box off,set(gca,'TickDir','out')
set(gca,'xtick',[])
set(gca,'FontSize',16)
ylabel('Power (Z)')
title('Miss Trials')
subplot(3,1,3),semilogy(plotWin,smoothdata(mean(LFP.narrowband.missLFP.betaGammaPvalue,1),'gaussian',1),'k','LineWidth',1),hold on
box off,set(gca,'TickDir','out')
set(gca,'FontSize',16)
xlabel('Time from cue onset (ms)') 
ylabel('P-value')
end
function output = itpcGEDstats(hit,miss)

hitprecue = max(hit(:,1:500),[],2);
hitpostcue = max(hit(:,500:1000),[],2);
hitpostRW = max(hit(:,1000:2000),[],2);
missprecue = max(miss(:,100:500),[],2);
misspostcue = max(miss(:,500:1000),[],2);
misspostRW = max(miss(:,1000:2000),[],2);
% plot
output.hitmprecue = mean(hitprecue);
output.hiteprecue = std(hitprecue)/sqrt(10);
output.hitmpostcue = mean(hitpostcue);
output.hitepostcue = std(hitpostcue)/sqrt(10);
output.hitmpostRW = mean(hitpostRW);
output.hitepostRW = std(hitpostRW)/sqrt(10);

output.missmprecue = mean(missprecue);
output.misseprecue = std(missprecue)/sqrt(10);
output.missmpostcue = mean(misspostcue);
output.missepostcue = std(misspostcue)/sqrt(10);
output.missmpostRW = mean(misspostRW);
output.missepostRW = std(misspostRW)/sqrt(10);


figure,errorbar(1:3,[output.hitmprecue,output.hitmpostcue,output.hitmpostRW],[output.hiteprecue,output.hitepostcue,output.hitepostRW],'b');hold on
errorbar(1:3,[output.missmprecue,output.missmpostcue,output.missmpostRW],[output.misseprecue,output.missepostcue,output.missepostRW],'k');
xlim([0 4]),ylim([0 3]),set(gca,'fontsize',16),box off
set(gca,'tickdir','out')
ylabel('Phase Synchonization (Z)')
legend('hit','miss')
end