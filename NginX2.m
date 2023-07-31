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
% make_UCLAChannelMap2(fpath,s); % Creates channel map for electrode array
make_UCLAChannelMap2(fpath)
kilosortPrep(Intan.allIntan,fpath)
set(0,'DefaultFigureWindowStyle','normal')
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
parameters.experiment = 'self'; % self - internally generated, cue - cue initiated
parameters.opto = 1; % 1 - opto ON , 0 - opto OFF
parameters.windowBeforePull = 1; % in seconds
parameters.windowAfterPull = 1; % in seconds
parameters.windowBeforeCue = 0.5; % in seconds
parameters.windowAfterCue = 1.5; % in seconds
parameters.Fs = 1000;
parameters.ts = 1/parameters.Fs;
%% Behavior
lfptime = 1/LFP.downSampleFreq:1/LFP.downSampleFreq:size(LFP.medianLFP,2)/LFP.downSampleFreq;
[Behaviour] = readLever3(parameters,lfptime);

%% Plot behaviour
for i=1:Behaviour.nHit
    plot(0:2000,smoothdata(Behaviour.hitTrace(i).trace),'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
    try
        hitTrace(i,:) = smoothdata(Behaviour.hitTrace(i).rawtrace);
    catch
        continue
    end
end
%%
buff = hitTrace;
buff(:,1:10) = 0;
buff (:,150:200) = 0;
 X = featureProject(smoothdata(buff,2,'gaussian',50),1,0,0);
%%
X = [];
count = 1;
for i = 1:31
    buff = smorate{i};
    blah = sum(buff,2);
    id = blah==0;
    buff(id,:) = [];
    normBuff = (buff-min(buff,[],'all'))/(max(buff,[],'all')-min(buff,[],'all'));
    if size(buff,1)>2
        X(:,:,count) = featureProject(normBuff,1,0,0);
        count = count+1;
    end
end
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
    hitWin = [Behaviour.hit(i,3)-1000, Behaviour.hit(i,3)+1000]; %intan position
    timestamphit(i,:) = [Behaviour.hit(i,4)-1, Behaviour.hit(i,4)+1]; %s
    hitLFP(:,:,i) = linearProbe(1:4,hitWin(1):hitWin(2));
    [waveletHit(:,:,i),fwavelet] = cwt(mean(hitLFP(:,:,i),1),1000,'FrequencyLimit',[4 80]);
    [powerCWThit(:,:,i), fwt] = calCWTSpectogram(mean(hitLFP(:,:,i),1),0:2000,1000,10,[4 80],0);
    [CSDoutputhit(:,:,i)]  = CSD(hitLFP(:,:,i)'/1E6,1000,20E-6);
end
for i = 1:Behaviour.nMiss
    missWin = [Behaviour.miss(i,3)-1000, Behaviour.miss(i,3)+1000]; %ms
    timestampmiss(i,:) = [Behaviour.miss(i,4)-1, Behaviour.miss(i,4)+1]; %s
    missLFP(:,:,i) = linearProbe(1:4,missWin(1):missWin(2));
    [waveletMiss(:,:,i),fwavelet] = cwt(mean(missLFP(:,:,i),1),1000,'FrequencyLimit',[4 80]);
    [powerCWTmiss(:,:,i), fwt] = calCWTSpectogram(mean(missLFP(:,:,i),1),0:2000,1000,10,[4 80],0);
    [CSDoutputmiss(:,:,i)]  = CSD(missLFP(:,:,i)'/1E6,1000,20E-6);
end
for i = 1:size(linearProbe,1)
tfmiss(:,:,i) = itpc(linearProbe(i,:),timestampmiss,1000,1);
tfhit(:,:,i) = itpc(linearProbe(i,:),timestamphit,1000,1);
end
%%
figure,plot(smoothdata(squeeze(mean(tfhit,1)),'gaussian',50),'k')
hold on,plot(smoothdata(squeeze(mean(tfmiss,1)),'gaussian',50),'r')

%%
[Shits,fhits,Serrhits]=mtspectrumc(squeeze(mean(hitLFP,1)),params);
[Smiss,fmiss,Serrmiss] = mtspectrumc(squeeze(mean(missLFP,1)),params);
%% Spikes analysis
path = [fpath,'/postAutoMerge'];
% Read in kilosort data for matlab analysis
SpikeClusters = readNPY(fullfile(path, 'spike_clusters.npy'));
SpikeSamples = readNPY(fullfile(path, 'spike_times.npy'));
Spikes.SpikeClusters = SpikeClusters; %Add one because of 0 index from python
Spikes.SpikeSamples = SpikeSamples;
Spikes = clusterSort(Spikes);
Spikes = ISI(Spikes,0.01,Fs,0); %Spikes, Interval, Fs
% Calculate Depth profile
load UCLA_chanMap
[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] =...
    spikeTemplatePosition(fpath,ycoords);
    for i = 1:length(tempAmps)
        Spikes.Clusters(i).spikeDepth = templateDepths(i);
        Spikes.Clusters(i).spikeAmplitude = tempAmps(i);
        Spikes.Clusters(i).waveforms = waveforms(i,:);
        Spikes.Clusters(i).spikeDuration = templateDuration(i)/Fs*1000;
    end
   %% Calculate PSTH
    count = 1;
    trials = {};
    for ii = 1:length(Spikes.Clusters)
        if ~isempty(Spikes.Clusters(ii).spikeTime)
            for i = 1:Behaviour.nHit
                trials{i,count}= zeros(1,length(0:0.001:Behaviour.hitTrace(i).t2-Behaviour.hitTrace(i).t1));
                spiketm = Spikes.Clusters(ii).spikeTime(Spikes.Clusters(ii).spikeTime>=Behaviour.hitTrace(i).t1 & Spikes.Clusters(ii).spikeTime<=Behaviour.hitTrace(i).t2);
                spiked = discretize(spiketm,Behaviour.hitTrace(i).t1:0.001:Behaviour.hitTrace(i).t2);
                if ~isnan(spiked)
                    trials{i,count}(spiked) = 1;
                end
                
            end
        else
            continue
        end
        count = count+1;
    end
    %% MISS
      %Calculate PSTH
    count = 1;
    trials = {};
    for ii = 1:length(Spikes.Clusters)
        if ~isempty(Spikes.Clusters(ii).spikeTime)
            for i = 1:Behaviour.nMiss
                trials{i,count}= zeros(1,length(0:0.001:Behaviour.missTrace(i).t2-Behaviour.missTrace(i).t1));
                spiketm = Spikes.Clusters(ii).spikeTime(Spikes.Clusters(ii).spikeTime>=Behaviour.missTrace(i).t1 & Spikes.Clusters(ii).spikeTime<=Behaviour.missTrace(i).t2);
                spiked = discretize(spiketm,Behaviour.missTrace(i).t1:0.001:Behaviour.missTrace(i).t2);
                if ~isnan(spiked)
                    trials{i,count}(spiked) = 1;
                end
                
            end
        else
            continue
        end
        count = count+1;
    end
    %%
    figure;
%****************
disp('Plotting rastergrams (slow) ...');
subplot('position',[0.1 0.4 0.4 0.55]);
make_nice_spike_raster(t);
V = axis;
%%
smorate = {}
t = {};
figure,hold on
grid on;
ylabel('Trial Number');
title(fprintf('Unit rasters'));
%*************
smooth_window = 100;  % give sigma of 12.5ms
for trils = 1:31
    count = 1;
    for i = 10:60
        t{1}= vertcat(trials{trils,i-5:i+5});
        try
        smorate{trils}(count,:) = make_nice_mean_raster(t,smooth_window);
        catch
            continue
        end
        count = count+1;
    end
end
V = axis;
%%
axis([0 size(spikea{1},2) V(3) V(4)]);
plot([interval(1),interval(1)],[V(3),V(4)],'k-'); hold on;
plot([interval(end),interval(end)],[V(3),V(4)],'k-'); hold on;
ylabel('Firing Rate');
xlabel('Time (ms)');
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
%% Some Functions
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