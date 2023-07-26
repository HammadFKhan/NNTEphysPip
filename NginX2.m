% clear; clc; 
% close all;
addpath(genpath('main'));
addpath(genpath('chronux'));
addpath(genpath('Kilosort'));
addpath(genpath('npy-matlab'));
addpath(genpath('spikes-master'));
IntanConcatenate
%%
% Intan = read_Intan_RHD2000_file(); %load intan data
% ParpoolConfig
fpath    = Intan.path; % where on disk do you want the analysis? ideally and SSD...
pathToYourConfigFile = strcat(pwd,'/main/'); % for this example it's ok to leave this path inside the repo, but for your own config file you *must* put it somewhere else!  
run(fullfile(pathToYourConfigFile, 'config_UCLAprobe.m'))
make_UCLAChannelMap2(fpath,s); % Creates channel map for electrode array
% make_UCLAMouseChannelMap(fpath);
kilosortPrep(Intan.allIntan,fpath)
set(0,'DefaultFigureWindowStyle','docked')
rez = KilosortAnalysis(fpath,ops);
% now fire up Phy and check these results. There should still be manual
% work to be done (mostly merges, some refinements of contaminated clusters). 
%% AUTO MERGES 
% after spending quite some time with Phy checking on the results and understanding the merge and split functions, 
% come back here and run Kilosort's automated merging strategy. This block
% will overwrite the previous results and python files. Load the results in
% Phy again: there should be no merges left to do (with the default simulation), but perhaps a few splits
% / cleanup. On realistic data (i.e. not this simulation) there will be drift also, which will usually
% mean there are merges left to do even after this step. 
% Kilosort's AUTO merges should not be confused with the "best" merges done inside the
% benchmark (those are using the real ground truth!!!)

%% LFP
set(0,'DefaultFigureWindowStyle','normal')
LFP = fastpreprocess_filtering(Intan.allIntan,10000);
LFP = bestLFP(LFP);
LFP = bandFilter(LFP,'depth'); % Extract LFPs based on 'depth' or 'single'
%% Parameters 
parameters.ephysFs = 1000;
parameters.ts = 1/parameters.ephysFs;
parameters.caTime = 0:parameters.ts:(size(LFP.medianLFP,2)-1)*parameters.ts;
parameters.windowBeforePull = 1; % in seconds
parameters.windowAfterPull = 1; % in seconds

%% Behavior
lfptime = 1/LFP.downSampleFreq:1/LFP.downSampleFreq:size(LFP.medianLFP,2)/LFP.downSampleFreq;
[Behaviour] = readLever(parameters,lfptime);
%% CSD and spectrogram
% Hit trials
CSDoutputhit = [];waveletHit = [];waveletMiss = [];powerCWThit = [];CSDoutputmiss = []; hitLFP = [];missLFP = [];
eIdx = find(s.sorted_probe_wiring(:,2)==0);
linearProbe = LFP.medianLFP(eIdx,:);
params.tapers = [5 9];
params.Fs = 1000;
params.fpass = [0 80];
params.err = [2 0.05];
for i = 1:100
    hitWin = [Behaviour.hit(i,3)-1000, Behaviour.hit(i,3)+1000]; %ms
    hitLFP(:,:,i) = linearProbe(1:4,hitWin(1):hitWin(2));
    [waveletHit(:,:,i),fwavelet] = cwt(mean(hitLFP(:,:,i),1),1000,'FrequencyLimit',[4 80]);
    [powerCWThit(:,:,i), fwt] = calCWTSpectogram(mean(hitLFP(:,:,i),1),0:2000,1000,10,[4 80],0);
    [CSDoutputhit(:,:,i)]  = CSD(hitLFP(:,:,i)'/1E6,1000,20E-6);
end
for i = 1:Behaviour.nMiss
    missWin = [Behaviour.hit(i,3)-1000, Behaviour.hit(i,3)+1000]; %ms
    missLFP(:,:,i) = linearProbe(1:4,missWin(1):missWin(2));
    [waveletMiss(:,:,i),fwavelet] = cwt(mean(missLFP(:,:,i),1),1000,'FrequencyLimit',[4 80]);
    [powerCWTmiss(:,:,i), fwt] = calCWTSpectogram(mean(missLFP(:,:,i),1),0:2000,1000,10,[4 80],0);
    [CSDoutputmiss(:,:,i)]  = CSD(missLFP(:,:,i)'/1E6,1000,20E-6);
end

[Shits,fhits,Serrhits]=mtspectrumc(squeeze(mean(hitLFP,1)),params);
[Smiss,fmiss,Serrmiss] = mtspectrumc(squeeze(mean(missLFP,1)),params);
%% Spikes
path = [fpath,'/postAutoMerge'];
% Read in kilosort data for matlab analysis
SpikeClusters = readNPY(fullfile(path, 'spike_clusters.npy'));
SpikeSamples = readNPY(fullfile(path, 'spike_times.npy'));
%% Analysis
Spikes.SpikeClusters = SpikeClusters+1; %Add one because of 0 index from python
Spikes.SpikeSamples = SpikeSamples;
Spikes = clusterSort(Spikes);
%% Calculate Depth profile
load UCLA_chanMap_fixed
[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] =...
    spikeTemplatePosition(fpath,ycoords);

%% Time-Frequency Analysis
[TimeFreq,LFP,betaGroup,Spikes] = tfAnalysis(Spikes,LFP,1); %Behavior state running 1 (0 rest)
[TimeFreq,LFP,betaGroupRest,Spikes] = tfAnalysis(Spikes,LFP,0,TimeFreq); %Behavior state running 1 (0 rest)
%% Save ITPC 
% Trunct rest and run into two columns for ease of comparison
L23thetaITPC = [TimeFreq.tfRest.depth.L23.theta.itpc TimeFreq.tfRun.depth.L23.theta.itpc];
L23betaITPC = [TimeFreq.tfRest.depth.L23.beta.itpc TimeFreq.tfRun.depth.L23.beta.itpc];
L23gammaITPC = [TimeFreq.tfRest.depth.L23.gamma.itpc TimeFreq.tfRun.depth.L23.gamma.itpc];

L5thetaITPC = [TimeFreq.tfRest.depth.L5.theta.itpc TimeFreq.tfRun.depth.L5.theta.itpc];
L5betaITPC = [TimeFreq.tfRest.depth.L5.beta.itpc TimeFreq.tfRun.depth.L5.beta.itpc];
L5gammaITPC = [TimeFreq.tfRest.depth.L5.gamma.itpc TimeFreq.tfRun.depth.L5.gamma.itpc];

%%
% [TimeFreq,LFP,betaGroupRest,Spikes] = tfAnalysis(Spikes,LFP,0,TimeFreq); %Behavior state running 1 (0 rest)

% plotTF(TimeFreq,LFP)
% TF stats of depth
TimeFreq.tf = TimeFreq.tfRun;
stats = tfStats(TimeFreq);ylim([0 0.4])
%%
tfDepth = TimeFreq.tf.depth;
betaGammaCoupling = gammaBetaCoupling(LFP,TimeFreq.tfRun,betaGroup);
betaGammaCouplingRest = gammaBetaCoupling(LFP,TimeFreq.tfRest,betaGroupRest);
betaGammam = mean(betaGammaCoupling,3);
figure,imagesc(-179:20:180,1:64,interp2(betaGammam')),colormap(jet)
figure,plot(mean(betaGammam))
% figure,imagesc(betaGammam(:,1:40)'),colormagithup(jet)
% figure,imagesc(betaGammam(:,41:64)'),colormap(jet)
%% Broadband Spectrograms per behavior state
pxSpecs = statePowerSpec(TimeFreq,LFP);
%% Plot Power Spectrum
f = pxSpecs.f;
broadband = pxSpecs.broadband;
mbroadband = mean(pxSpecs.broadband,2);
figure, hold on,semilogx(pxSpecs.f(2:end),pxSpecs.broadband(2:end,:)),semilogx(pxSpecs.f,mbroadband,'k','LineWidth',3),xlim([0 100]),set(gca, 'XScale', 'log')
figure,lineError(pxSpecs.f(2:end),pxSpecs.broadband(2:end,:)','ste'),set(gca, 'XScale', 'log'),xlim([0 100]),ylim([-100 40]),box off
figure,imagesc(-1000:1000,pxSpecs.fwav,pxSpecs.waveletNorm),colormap(jet),axis xy,colorbar,caxis([-.65 .65])
%Segment theta,beta,gamma from power spectrum
thetaP = pxSpecs.broadband((pxSpecs.f>4 & pxSpecs.f<12),:);
betaP = pxSpecs.broadband((pxSpecs.f>12 & pxSpecs.f<32),:);
gammaP1 = pxSpecs.broadband((pxSpecs.f>32 & pxSpecs.f<58),:); % Cutout 60 Hz
gammaP2 =  pxSpecs.broadband((pxSpecs.f>62 & pxSpecs.f<100),:);
gammaP = [gammaP1;gammaP2];
idx1 = mean(thetaP,1);idx2 = mean(betaP,1);idx3 = mean(gammaP,1);
total = [mean(idx1) mean(idx2) mean(idx3)];
err = [std(idx1)/sqrt(length(idx1)) std(idx2)/sqrt(length(idx2)) std(idx3)/sqrt(length(idx3))];
figure,bar(total),hold on
errorbar(1:3,total,err)
idx1 = idx1';idx2 = idx2';idx3 = idx3';
%% 
[L23RunAvg,L5RunAvg,L23UnitRestFreq,L5UnitRestFreq] = spikeRaster(Spikes,1,5); %Spikes, flag (0/1) for rest/run, and scatter size (sz)
[L23RestAvg,L5RestAvg,L23UnitRestFreq,L5UnitRestFreq] = spikeRaster(Spikes,0,5); %Spikes, flag (0/1) for rest/run
%%
t1 = cellfun('size',L23RestAvg,1)/4;
t2 = cellfun('size',L5RestAvg,1)/4;
L23Rest = t1';L5Rest = t2';
figure,histogram(t1,0:10:100),hold on
histogram(t2,0:10:100), box off

t = vertcat(t1',t2');
x = [repmat({'L23'},length(t1),1);];
x1 = repmat({'L5'},length(t2),1);
figure,boxplot(t1,x,'plotstyle','compact'), hold on, box off,ylim([0 100]),title('L23 Rest')
figure,boxplot(t2,x,'plotstyle','compact'), hold on, box off,ylim([0 100]),title('L5 Rest')


t1 = cellfun('size',L23RunAvg,1)/4;
t2 = cellfun('size',L5RunAvg,1)/4;
L23Run = t1';L5Run = t2';
Resttotal = [L23Rest,L5Rest];
Runtotal = [L23Run,L5Run];
%%
figure,histogram(t1,0:10:150),hold on
histogram(t2,0:10:150), box off

t = vertcat(t1',t2');
%x = [repmat({'L23'},length(t1),1);repmat({'L5'},length(t2),1);];
x = [repmat({'L23'},length(t1),1);];
x1 = repmat({'L5'},length(t2),1);
figure,boxplot(t1,x,'plotstyle','compact'), hold on, box off,ylim([0 100]),title('L23 Run')
figure,boxplot(t2,x1,'plotstyle','compact'), hold on, box off,ylim([0 100]),title('L5 Run')

%% Beta Analysis for each electrode
if exist('bstats','var')
    clear bstats
    clear betaNorm
    clear csd
    clear peakAlign
end
% badElectrode = 44;
electrode = 1:size(betaGroup,2);
% electrode(badElectrode) = [];
for i = electrode % Checks electrode size for median
    disp(['Electrode: ' num2str(i)])
    [peakAlign{i},mLFP{i},betaNorm{i},bstats(i)] = betaAnalysis(betaGroup(i).electrode,LFP.LFP);
end 


% Beta Event Rate
betaEventRate = [];
betaEventDuration = [];
for i = electrode
    betaEventRate = [betaEventRate; bstats(i).betaER];
    betaEventDuration= [betaEventDuration;bstats(i).betaDuration];
end
figure,boxplot(betaEventRate,'Plotstyle','compact'),ylim([1 6]),box off
figure,boxplot(betaEventDuration,'Plotstyle','compact'),box off
figure,histogram(betaEventRate,1:15), box off, set(gca,'TickDir','out');
figure,histogram(betaEventDuration,25:10:250),box off, set(gca,'TickDir','out');
%%
% Take out non-existant cell fields
betaNorm = betaNorm(~cellfun('isempty',betaNorm));
peakAlign = peakAlign(~cellfun('isempty',peakAlign));
csd = csd(~cellfun('isempty',csd));
%% Beta stats
set(0,'DefaultFigureWindowStyle','normal')
norm = vertcat(betaNorm{:});

for i = 1:length(betaNorm)
    peakNorm(i,:) = max(betaNorm{i},[],2);
end
figure,imagesc(f,1:64,interp2(peakNorm)),colormap(jet),axis xy,set(gcf, 'Position',  [100, 100,300, 500])
figure,imagesc(-100:100,-100:100,norm),colormap(jet),colorbar,axis xy,set(gcf, 'Position',  [100, 100, 300, 500])
% Plot stats across electrodes
LFPdepth = linspace(1,round(Spikes.Depth.depth(end),-2),64); %Round to nearest 100th micron
if isstruct(bstats)
    bstatsOrg = bstats;
    bstats = (struct2cell(bstats));
    bstats = squeeze(bstats)';
end
% figure('Name','Beta Duration'),bar(LFPdepth,bstats(:,1),'BarWidth',1),set(gcf, 'Position',  [100, 100, 500, 500])
% figure('Name','Beta Amplitude'),bar(LFPdepth,bstats(:,2),'BarWidth',1),set(gcf, 'Position',  [100, 100, 500, 500])
% figure('Name','Beta Event Rate'),bar(LFPdepth,bstats(:,3),'BarWidth',1),set(gcf, 'Position',  [100, 100, 500, 500])
% Bar plot for each layer
stats = betaStats(bstats,LFPdepth,0); %(Bstats, LFPdepth,plotFlag

%% Plot beta traces for each electrode
% figure, hold on
for i = 1:size(peakAlign,2) % Checks electrode size for median
    mPeakAlign(:,i) = mean(peakAlign{i},1);
end
mPeakAlign = mPeakAlign(:,[1:43,45:64]);

figure,plot(mean(mPeakAlign,2)), hold on
plot(mean(mPeakAlign,2)+std(mPeakAlign,0,2),'--r')
plot(mean(mPeakAlign,2)-std(mPeakAlign,0,2),'--r')

figure,stack_plot(flip(mPeakAlign'),0.35,1)
normPeakAlign = (mPeakAlign-min(mPeakAlign,[],'all'))/(max(mPeakAlign,[],'all')-min(mPeakAlign,[],'all'));
figure,imagesc(0:250,LFPdepth,interp2(smoothdata(normPeakAlign'))),caxis([.1 1])
%% Plot beta CSD for each electrode
for i = 1:size(csd,2) % Checks electrode size for median
    mcsd(:,:,i) = mean(csd{i},3);
end
mcsd = mean(mcsd,3);
figure,imagesc(0:250,LFPdepth,interp2(smoothdata((mcsd')),2)),colormap(jet),caxis([-0.2 0.2])

%%
behaviorflag = 1;
spikeTriggeredBeta = betaEventPSH(betaGroup,Spikes,behaviorflag); %set behavior flag 0 or 1 for rest/run
