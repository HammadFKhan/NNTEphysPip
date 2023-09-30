% clear; clc; 
% close all;
addpath(genpath('Main'));
addpath(genpath('chronux'));
% addpath(genpath('Kilosort'));
addpath(genpath('npy-matlab'));
addpath(genpath('spikes-master'));
% IntanConcatenate legacy version
ds_filename = intanPreprocessing2; %% double check file type
% Run Kilosort3 
% load only neccessary variables from memory mapped file
data = matfile(ds_filename);
fpath = data.fpath;
% Kilosort264FTestcode
savepath = fullfile(fpath,['loadme','.mat']);
save(savepath,'ds_filename');
clearvars -except ds_filename
%% Parameters for behaviour
data = matfile(ds_filename);
parameters.experiment = 'cue'; % self - internally generated, cue - cue initiated
parameters.opto = 0; % 1 - opto ON , 0 - opto OFF
parameters.windowBeforePull = 1; % in seconds
parameters.windowAfterPull = 1; % in seconds
parameters.windowBeforeCue = 0.5; % in seconds
parameters.windowAfterCue = 1.5; % in seconds
parameters.Fs = 1000; % Eventual downsampled data
parameters.ts = 1/parameters.Fs;
[Behaviour] = readLever(parameters,data.amplifierTime);
if size(data.digitalChannels,1)>3
    [IntanBehaviour] = readLeverIntan(parameters,data.amplifierTime,data.analogChannels,data.digitalChannels(2:end,:),Behaviour);
else
    [IntanBehaviour] = readLeverIntan(parameters,data.amplifierTime,data.analogChannels,data.digitalChannels,Behaviour);
end
% Calculate ITI time for trials and reward/no reward sequence
temp1 = arrayfun(@(x) x.LFPtime(1), IntanBehaviour.cueHitTrace);
temp1 = vertcat(temp1,ones(1,IntanBehaviour.nCueHit)); %  write 1 for reward given
temp2 = arrayfun(@(x) x.LFPtime(1), IntanBehaviour.cueMissTrace);
temp2 = vertcat(temp2,zeros(1,IntanBehaviour.nCueMiss)); %  write 0 for no reward given
temp = [temp1,temp2];
[~,idx] = sort(temp(1,:)); %sort by occurance
IntanBehaviour.ITI = temp(:,idx);
%% Plot behaviour
figure
for i=1:IntanBehaviour.nCueHit
    plot(0:2000,smoothdata(IntanBehaviour.cueHitTrace(i).trace),'Color',[0 0 0 0.2],'LineWidth',1.5);
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
load UCLA_chanMap_64F2
% load UCLA_chanmap_fixed.mat
if ~exist('lfp','var'),lfp = data.amplifierData;end
%TODO check if the field orientation during insertion is reversed (ie. probe 1 is lateral to probe 2)
% probe1 = lfp(s.sorted_probe_wiring(:,5)==1,:);
% probe2 = lfp(s.sorted_probe_wiring(:,5)==2,:);
% chanProbe1 = s.sorted_probe_wiring(s.sorted_probe_wiring(:,5)==1,:); %needed for linear channel mapping later
% chanProbe2 = s.sorted_probe_wiring(s.sorted_probe_wiring(:,5)==2,:); 
% clear lfp
probe1 = lfp;
% LFP filter
set(0,'DefaultFigureWindowStyle','normal')
LFP.probe1= fastpreprocess_filtering(probe1,data.targetedFs);
LFP.probe1 = bestLFP(LFP.probe1);
% LFP.probe1 = bandFilter(LFP.probe1,'depth'); % Extract LFPs based on 'depth' or 'single'
% LFP.probe2 = fastpreprocess_filtering(probe2,data.targetedFs);
% LFP.probe2 = bestLFP(LFP.probe2);
% % LFP.probe2 = bandFilter(LFP.probe2,'depth'); % Extract LFPs based on 'depth' or 'single'
% % Build linear channels (should be four from the 64F)
% LFP.probe1.chan1 = find(chanProbe1(:,2)>10);
% LFP.probe1.chan2 = find(chanProbe1(:,2)==-10);
% LFP.probe2.chan1 = find(chanProbe2(:,2)==290.1);
% LFP.probe2.chan2 = find(chanProbe2(:,2)==310.1);
%% CSD and spectrogram
LFP.probe1.spectrogram = leverLFPAnalysis(LFP.probe1.LFP,IntanBehaviour);
% LFP.probe2.spectrogram = leverLFPAnalysis(LFP.probe2.LFP,IntanBehaviour);
%% Now we used the segmented LFPs to do statistical analysis across narrowbands
LFP.probe1.narrowband.hitLFP = bandFilter2(LFP.probe1.spectrogram.hitLFP,1000);
% LFP.probe2.narrowband.hitLFP = bandFilter2(LFP.probe2.spectrogram.hitLFP,1000);
LFP.probe1.narrowband.missLFP = bandFilter2(LFP.probe1.spectrogram.missLFP,1000);
% LFP.probe2.narrowband.missLFP = bandFilter2(LFP.probe2.spectrogram.missLFP,1000);
%% Calculate trial-trial narrow power using hilbert across electrode
for electrode = 1:32
    % Beta
    [LFP.probe1.narrowband.hitLFP.betaPower(:,:,electrode),LFP.probe1.narrowband.missLFP.betaPower(:,:,electrode)]=...
        calcBandPower(squeeze(LFP.probe1.narrowband.hitLFP.beta_band(electrode,:,:)),squeeze(LFP.probe1.narrowband.missLFP.beta_band(electrode,:,:)));
    
    % Gamma
    [LFP.probe1.narrowband.hitLFP.gammaPower(:,:,electrode),LFP.probe1.narrowband.missLFP.gammaPower(:,:,electrode)] =...
        calcBandPower(squeeze(LFP.probe1.narrowband.hitLFP.gamma_band(electrode,:,:)),squeeze(LFP.probe1.narrowband.missLFP.gamma_band(electrode,:,:)));
end
%% Calculate p-value of beta/gamma coupling for each electrode pair
for nn = 1:size(LFP.probe1.narrowband.hitLFP.betaPower,3) %electrode
    for n = 1:size(LFP.probe1.narrowband.hitLFP.betaPower,2) %time
        LFP.probe1.narrowband.hitLFP.betaGammaPvalue(nn,n) = ranksum(LFP.probe1.narrowband.hitLFP.betaPower(:,n,nn),LFP.probe1.narrowband.hitLFP.gammaPower(:,n,nn));
        LFP.probe1.narrowband.missLFP.betaGammaPvalue(nn,n) = ranksum(LFP.probe1.narrowband.missLFP.betaPower(:,n,nn),LFP.probe1.narrowband.missLFP.gammaPower(:,n,nn));
    end
end
%% Beta vs Gamma Power
plotBetaGammaPower(IntanBehaviour.parameters,LFP.probe1) %Reference probe here
%% Plot Power by depth
sz = size(LFP.probe1.narrowband.hitLFP.gammaPower);
temp = smoothdata((squeeze(mean(LFP.probe1.narrowband.hitLFP.betaPower,1))),'gaussian',1);
figure,
imagesc(-IntanBehaviour.parameters.windowBeforeCue*1000:IntanBehaviour.parameters.windowAfterCue*1000,1:sz(3),temp'),colormap(hot)
xlabel('Time from cue (ms)')
ylabel('Depth')
set(gca,'fontsize',16)
%% Beta vs gamma phase synchonization p-value
idxBeta = find(LFP.probe1.spectrogram.frex>12 & LFP.probe1.spectrogram.frex<30); %Frequency regime
idxGamma = find(LFP.probe1.spectrogram.frex>30 & LFP.probe1.spectrogram.frex<80); %Frequency regime
for nn = 1:size(LFP.probe1.spectrogram.tfhit,3) %electrode
    for n = 1:size(LFP.probe1.spectrogram.tfhit,2) %timepoint
        LFP.probe1.narrowband.hitLFP.phasebetaGammaPvalue(nn,n) = ...
            ranksum(squeeze(LFP.probe1.spectrogram.tfhit(idxBeta,n,nn)),LFP.probe1.spectrogram.tfhit(idxGamma,n,nn));
        LFP.probe1.narrowband.missLFP.phasebetaGammaPvalue(nn,n) = ...
            ranksum(LFP.probe1.spectrogram.tfmiss(idxBeta,n,nn),LFP.probe1.spectrogram.tfmiss(idxGamma,n,nn));
    end
end
%% plot ITPC across depth
broadBandIdx = 7:62; %Broadband idx base on spectrogram.frex
sz = size(LFP.probe1.spectrogram.tfhit);
temp = smoothdata(squeeze(mean(LFP.probe1.spectrogram.tfhit(broadBandIdx,:,:),1)),'gaussian',10);
temp = smoothdata(temp,2);
figure,
imagesc(-IntanBehaviour.parameters.windowBeforeCue*1000:IntanBehaviour.parameters.windowAfterCue*1000,1:sz(3),temp'),colormap(hot)
xlabel('Time from cue (ms)')
ylabel('Depth')
set(gca,'fontsize',16)
%% Prestim and post stim beta/gamma
% Here we plot some sort of correlation map between precue response
% dynamics and underlying reaction time for lever pull
% Beta hit
[LFP.probe1.narrowband.hitLFP.preBetastim,LFP.probe1.narrowband.hitLFP.postBetastim,LFP.probe1.narrowband.hitLFP.preBetastimIdx,LFP.probe1.narrowband.hitLFP.postBetastimIdx] =...
    stimulusPower(LFP.probe1.narrowband.hitLFP.betaPower,IntanBehaviour);
%Gamma hit
[LFP.probe1.narrowband.hitLFP.preGammastim,LFP.probe1.narrowband.hitLFP.postGammastim,LFP.probe1.narrowband.hitLFP.preGammastimIdx,LFP.probe1.narrowband.hitLFP.postGammastimIdx] =...
    stimulusPower(LFP.probe1.narrowband.hitLFP.gammaPower,IntanBehaviour);
% Beta miss
[LFP.probe1.narrowband.missLFP.preBetastim,LFP.probe1.narrowband.missLFP.postBetastim,LFP.probe1.narrowband.missLFP.preBetastimIdx,LFP.probe1.narrowband.missLFP.postBetastimIdx] =...
    stimulusPower(LFP.probe1.narrowband.missLFP.betaPower,IntanBehaviour);
% Gamma miss
[LFP.probe1.narrowband.missLFP.preGammastim,LFP.probe1.narrowband.missLFP.postGammastim,LFP.probe1.narrowband.missLFP.preGammastimIdx,LFP.probe1.narrowband.missLFP.postGammastimIdx] =...
    stimulusPower(LFP.probe1.narrowband.missLFP.gammaPower,IntanBehaviour);
%% Phase amplitude coupling prestim and poststim/ hit and miss
phas_freqs = 8:2:30;
ampl_freqs = 30:4:80;
count = 1;
preCueHitphaseamp = [];
postCueHitphaseamp = [];
preCueMissphaseamp = [];
postCueMissphaseamp = [];
electrode = 32;
for n = 1:196 %trial
    preCueHitphaseamp(:,:,count) = phaseAmpCoupling(squeeze(compts(1,1:500,n)),1000,phas_freqs,ampl_freqs);
    postCueHitphaseamp(:,:,count) = phaseAmpCoupling(squeeze(compts(1,501:2001,n)),1000,phas_freqs,ampl_freqs);
%     preCueMissphaseamp(:,:,count) = phaseAmpCoupling(LFP.probe1.spectrogram.missLFP(electrode,1:500,n),1000,phas_freqs,ampl_freqs);
%     postCueMissphaseamp(:,:,count) = phaseAmpCoupling(LFP.probe1.spectrogram.missLFP(electrode,501:2000,n),1000,phas_freqs,ampl_freqs);
    count = count+1;
end


%%
% z-score normalized to shuffled
n_iter = 200;
% bm = zeros(1,length(n_iter));
for bi=1:n_iter
    cutpoint = randsample(round(size(hitCuephaseamp,3)/10):round(size(hitCuephaseamp,3)*.9),6);
    bm(:,:,bi) = mean(hitCuephaseamp(:,:,cutpoint),3);
end
temp = mean(hitCuephaseamp,3)-mean(bm,3)./std(bm,[],3);

for bi=1:n_iter
    cutpoint = randsample(round(size(missCuephaseamp,3)/10):round(size(missCuephaseamp,3)*.9),6);
    bm(:,:,bi) = mean(missCuephaseamp(:,:,cutpoint),3);
end
temp = mean(missCuephaseamp,3)-mean(bm,3)./std(bm,[],3);


%% Plot
figure(1),clf,contourf(phas_freqs,ampl_freqs,preCueHitphaseamp(:,:,1)',140,'linecolor','none'),colormap(jet)
set(gca,'clim',[-2.98 2.98]),colormap(jet)
xlabel('Frequency for phase')
ylabel('frequency for amplitude')
title('Pre-Cue Hit phase-amplitude coupling')
figure(2),clf,contourf(phas_freqs,ampl_freqs,mean(postCueHitphaseamp,3)',140,'linecolor','none'),colormap(jet)
set(gca,'clim',[-.5 .5]),colormap(jet)
xlabel('Frequency for phase')
ylabel('frequency for amplitude')
title('Post-cue Hit phase-amplitude coupling')
figure(3),clf,contourf(phas_freqs,ampl_freqs,mean(preCueMissphaseamp,3)',140,'linecolor','none'),colormap(jet)
set(gca,'clim',[-.5 .5]),colormap(jet)
xlabel('Frequency for phase')
ylabel('frequency for amplitude')
title('Pre-cue Miss phase-amplitude coupling')
figure(4),clf,contourf(phas_freqs,ampl_freqs,mean(postCueMissphaseamp,3)',140,'linecolor','none'),colormap(jet)
set(gca,'clim',[-.5 .5]),colormap(jet)
xlabel('Frequency for phase')
ylabel('frequency for amplitude')
title('Post-cue Miss phase-amplitude coupling')
%% gedCFA
% Calculate channel coherence
rawData = LFP.probe1.LFP(:,:);
Rdata = LFP.probe1.spectrogram.hitLFP(:,1:1500,:);
Sdata = LFP.probe1.spectrogram.hitLFP(:,1501:3000,:);
GED = genEigenDecomp(rawData,Rdata,Sdata);
%% Calculate generalize phase for compts
xo = bandpass_filter(GED.compts,5,40,1000); %x,f1,f2,Fs
sz = size(xo);
xo = reshape(xo,sz(1),1,sz(2));
[xgp,wt] = generalized_phase(xo, 1000, 0 );
%%
for i = 1:IntanBehaviour.nCueHit
    disp(['Analyzing trial: ' num2str(i)])
    hitWin = floor(IntanBehaviour.cueHitTrace(i).LFPtime*1000); %multiply by Fs;
    hitxgp{i} = xgp(:,hitWin);
end
for i = 1:80
    disp(['Analyzing trial: ' num2str(i)])
    missWin = floor(IntanBehaviour.cueMissTrace(i).LFPtime*1000);
    missxgp{i} = xgp(:,missWin);
end

%%
hitxgp = cellfun(@(x) reshape(x,3,1,[]),hitxgp,'UniformOutput',false);
missxgp = cellfun(@(x) reshape(x,3,1,[]),missxgp,'UniformOutput',false);
%%
[PA] = calPhaseAlignment(missxgp);
figure,plot((squeeze(PA)'),'LineWidth',2);
cmap = (gray(4));
set(gca(),'ColorOrder',cmap)
%% Now do the same analysis of LFP but using broadband GED LFP
LFP.probe1.GED.comp1 = leverLFPAnalysis(GED.compts(1,:),IntanBehaviour);
LFP.probe1.GED.comp2 = leverLFPAnalysis(GED.compts(2,:),IntanBehaviour);
LFP.probe1.GED.comp3 = leverLFPAnalysis(GED.compts(3,:),IntanBehaviour);
%% ITPC by component
LFP.probe1.GED.itpc = GEDITPC(LFP.probe1.GED);
figure,
subplot(131),plotSpectrogram(LFP.probe1.GED.comp1.tfhit,-1499:1500,LFP.probe1.GED.comp1.frex,'contourf'); caxis([-3 3])
subplot(132),plotSpectrogram(LFP.probe1.GED.comp2.tfhit,-1499:1500,LFP.probe1.GED.comp1.frex,'contourf'); caxis([-3 3])
subplot(133),plotSpectrogram(LFP.probe1.GED.comp3.tfhit,-1499:1500,LFP.probe1.GED.comp1.frex,'contourf'); caxis([-3 3])
figure,subplot(131),plot(-1499:1500,LFP.probe1.GED.itpc.comp1.hit)
subplot(132),plot(-1499:1500,LFP.probe1.GED.itpc.comp2.hit)
subplot(133),plot(-1499:1500,LFP.probe1.GED.itpc.comp3.hit)
% %%% Statistically pool 
% [val,idx] = max(LFP.probe1.GED.itpc.comp1.hit);
% idx = idx-500;
%%
hit = LFP.probe1.GED.itpc.comp3.hit;
miss = LFP.probe1.GED.itpc.comp3.miss;
figure,hold on
plot(-1499:1500,hit,'k','LineWidth',2)
plot(-1499:1500,miss,'Color',[188/255 190/255 192/255],'LineWidth',2)
xlabel('Time from cue (ms)')
set(gca,'fontsize',16)
box off, set(gca,'TickDir','out')
yline(1.98)
xline(mean(IntanBehaviour.reactionTime)*1000)
xline(0)
xlim([-500 1500])
%% Analyze as a function of phase locking
%
freq4phase = 15; % in Hz
freq4power = 40; 
for n = 1:100
% get phase values
phasefilt = filterFGx(squeeze(LFP.probe1.GED.comp1.hitLFP(:,:,n)),1000,freq4phase,2);
phase(n,:) = angle(hilbert(phasefilt));

% get power values (note: 'power' is a built-in function so we'll name this variable 'amp')
ampfilt = filterFGx(squeeze(LFP.probe1.GED.comp1.hitLFP(:,:,1)),1000,freq4power,40,0);
amp(n,:) = abs(hilbert(ampfilt)).^2;

% % plot power and phase
% figure(2), clf
% plot(phase)
% hold on
% plot((amp-mean(amp))/std(amp),'r')
% hint: zoom in on part of these data!

%plot power as a function of phase
figure(3), 
polarplot(mean(phase,1),mean(amp,1),'.')
n_hist_bins = 30;

% phase_edges = linspace(min(phase),max(phase),n_hist_bins+1);
% amp_by_phases = zeros(1,n_hist_bins);
% 
% for i=1:n_hist_bins-1
%     amp_by_phases(i) = mean(amp(phase>phase_edges(i) & phase<phase_edges(i+1)));
% end

% figure(4),hold on
% bar(phase_edges(1:end-1),amp_by_phases,'histc');
% set(gca,'xlim',[-pi-.5 pi+.5])
% xlabel([ 'Phase at ' num2str(freq4phase) ' Hz (rad.)' ])
% ylabel([ 'Power at ' num2str(freq4power) ' Hz' ])
% box off,set(gca,'TickDir','out')
% set(gca,'fontsize',16)
end
% significance through permutation testing 

%% compute amplitude-modulated phase
observed_modulation = abs(mean(amp.*exp(1i*phase)));

% now run permutation testing to see the significance of this interaction
num_iterations = 10000;
permuted_results = zeros(1,num_iterations);
npnts = length(amp);

for ni=1:num_iterations
    % randomly resort power values
    cutpoint = randsample(round(npnts/10):round(npnts*.9),1);
    permuted_results(ni) = abs(mean(amp([ cutpoint:end 1:cutpoint-1 ]).*exp(1i*phase)));
end

% print value
z_modulation = (observed_modulation-mean(permuted_results))/std(permuted_results);
disp([ 'Z-modulation of ' num2str(freq4power) ' Hz power modulated by ' num2str(freq4phase) ' Hz phase is ' num2str(z_modulation) '.' ])

% show histogram of permuted and real values
figure(5), clf
histogram(permuted_results,50);
hold on
plot([observed_modulation observed_modulation],get(gca,'ylim')/2,'m-p','linewi',4,'markersize',16);
legend({'histogram of permuted values';'observed value'})
xlabel('Absolute modulate (arb. units depending on power values)')
ylabel('Count')


%% Spikes analysis
data = matfile(ds_filename);
path = [data.fpath,'/kilosort3/merged'];
%%% Spike preprocessing (includes merging (optional) and channel info
%%% return)
% Read in kilosort data for matlab analysis
SpikeClusters = readNPY(fullfile(path, 'spike_clusters.npy'));
SpikeSamples = readNPY(fullfile(path, 'spike_times.npy'));
SpikeChannel = readNPY(fullfile(path,'channel_positions.npy'));
% Kilosort3AutoMergeTester
Spikes.SpikeClusters = SpikeClusters; 
Spikes.SpikeSamples = SpikeSamples;
Spikes = clusterSort(Spikes);
Spikes = ISI(Spikes,0.01,data.Fs,0); %Spikes, Interval, Fs
% Calculate Depth profile
load chanMap64F2
% load UCLA_chanMap
[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms, max_site] =...
    spikeTemplatePosition(data.fpath,ycoords,'invert');
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
%% Basic spike analysis
% z-score spike rates
[Spikes] = sortSpkLever(Spikes,IntanBehaviour);
%% Make it layer specific TODO:Gamma GED/Spike coherence
Spikes = layerspikeSort(Spikes);
%%
figure,hold on
for n = 61
subplot(2,1,[1]),Show_Spikes(Spikes.PSTH.hit.spks{n}),axis off
subplot(2,1,[2]),bar(-500:1500,smoothdata(Spikes.PSTH.hit.spkRates(n,:)),'FaceColor',[28/255 117/255 188/255],'EdgeColor','none')
axis tight, box off, set(gca,'TickDir','out')
set(gca,'fontsize',16)
xline(0)
xline(mean(IntanBehaviour.reactionTime)*1000)
end
%%
figure,hold on
for n = 3
subplot(2,1,[1]),Show_Spikes(Spikes.PSTH.hit{n}),axis off
subplot(2,1,[2]),bar(-1500:500,smoothdata(Spikes.PSTH.hitspikeRates(n,:)),'FaceColor',[169/255 79/255 25/255],'EdgeColor','none')
axis tight, box off, set(gca,'TickDir','out')
set(gca,'fontsize',16')
xline(0)
xline(mean(IntanBehaviour.reactionTime)*1000)
end
%%
figure,hold on
for n = 4
subplot(2,1,[1]),Show_Spikes(Spikes.PSTH.hit{n}),axis off
subplot(2,1,[2]),bar(-500:1500,smoothdata(Spikes.PSTH.hitspikeRates(n,:)),'FaceColor',[190/255 30/255 45/255],'EdgeColor','none')
axis tight, box off, set(gca,'TickDir','out')
set(gca,'fontsize',16')
xline(0)
xline(mean(IntanBehaviour.reactionTime)*1000)
end
%%
figure,hold on
for n = 33
subplot(2,1,[1]),Show_Spikes(Spikes.PSTH.hit{n}),axis off
subplot(2,1,[2]),bar(-500:1500,smoothdata(Spikes.PSTH.hitspikeRates(n,:)),'FaceColor',[43/255 57/255 144/255],'EdgeColor','none')
axis tight, box off, set(gca,'TickDir','out')
set(gca,'fontsize',16')
xline(0)
xline(mean(IntanBehaviour.reactionTime)*1000)
end
%% Neural Trajectory Analysis using GPFA
Spikes = makeSpikeGPFA(Spikes);
[result,seqTrain] = gpfaAnalysis(Spikes.GPFA.hit.dat,0); %Run index

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
%     [CSDoutputhit(:,:,i)]  = CSD(hitLFP(:,:,i)/1E6,1000,50E-6);
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
    %     [CSDoutputmiss(:,:,i)]  = CSD(missLFP(:,:,i)'/1E6,1000,20E-6);
end
for i = 1:size(linearProbe,1)
    disp(['Analyzing electrode: ' num2str(i)])
    tfmiss(:,:,i) = itpc(linearProbe(i,:),timestampmiss,1000,0);
    [tfhit(:,:,i),frex,pnts] = itpc(linearProbe(i,:),timestamphit,1000,0);
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
broadBandIdx = 30:54; %Broadband idx base on spectrogram.frex (35:40),30:35 low Beta, 55:60 gamma
temp = GED.comp1;
%Search freq space and note maximal response for beta
betatfhit = temp.tfhit(broadBandIdx,:);
[~,idx] = max(betatfhit,[],1);
idx = round(mean(idx));
tFreq = broadBandIdx(idx);
output.comp1.hit = smoothdata(squeeze(mean(temp.tfhit(tFreq-1:tFreq+1,:),1)));
output.comp1.miss = smoothdata(squeeze(mean(temp.tfmiss(tFreq-1:tFreq+1,:),1)));

broadBandIdx = 30:54; %Broadband idx base on spectrogram.frex (35:40),30:35 low Beta, 55:60 gamma
temp = GED.comp2;
betatfhit = temp.tfhit(broadBandIdx,:);
[~,idx] = max(betatfhit,[],1);
idx = round(mean(idx));
tFreq = broadBandIdx(idx);
output.comp2.hit = smoothdata(squeeze(mean(temp.tfhit(tFreq-1:tFreq+1,:),1)));
output.comp2.miss = smoothdata(squeeze(mean(temp.tfmiss(tFreq-1:tFreq+1,:),1)));

broadBandIdx = 30:54; %Broadband idx base on spectrogram.frex (35:40),30:35 low Beta, 55:60 gamma
temp = GED.comp3;
betatfhit = temp.tfhit(broadBandIdx,:);
[~,idx] = max(betatfhit,[],1);
idx = round(mean(idx));
tFreq = broadBandIdx(idx);
output.comp3.hit = smoothdata(squeeze(mean(temp.tfhit(tFreq-1:tFreq+1,:),1)));
output.comp3.miss = smoothdata(squeeze(mean(temp.tfmiss(tFreq-1:tFreq+1,:),1)));
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

function xo = bandpass_filter(x,f1,f2,Fs)
ct = [f1 f2];
ct = ct / (Fs/2);
filter_order = 4;
[b,a] = butter( filter_order, ct );
% proceed with filtering
for rr = 1:size(x,1) % electrode
        if ~all( isnan(squeeze(x(rr,:))) )
            xo(rr,:) = filtfilt( b, a, squeeze(x(rr,:)) );
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