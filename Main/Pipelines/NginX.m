% clear; clc; 
% close all;
addpath(genpath('main'));
addpath(genpath('chronux'));
addpath(genpath('Kilosort'));
addpath(genpath('npy-matlab'));
addpath(genpath('spikes-master'));
IntanConcatenate
% Intan = read_Intan_RHD2000_file(); %load intan data
global useGPU 
useGPU = 0;
% ParpoolConfig
fpath    = Intan.path; % where on disk do you want the analysis? ideally and SSD...
pathToYourConfigFile = strcat(pwd,'/main/'); % for this example it's ok to leave this path inside the repo, but for your own config file you *must* put it somewhere else!  
run(fullfile(pathToYourConfigFile, 'config_eMouse.m'))
%make_UCLAChannelMap2(fpath,s); % Creates channel map for electrode array
make_UCLAMouseChannelMap(fpath);
%%
% filtData = preprocess_filtering(Intan.allIntan(:,1:400000),Intan.t_amplifier);
%% Ripples
% [Ripples,filtData] = SWR(Intan);
% plotClusterless(Ripples,filtData,Intan)


%% Kilosort Analysis
useGPU = 1;
kilosortPrep(Intan.allIntan,fpath)
set(0,'DefaultFigureWindowStyle','normal')
rez = KilosortAnalysis(fpath,ops);
% now fire up Phy and check these results. There should still be manual
% work to be done (mostly merges, some refinements of contaminated clusters). 
useGPU = 0;

%% AUTO MERGES 
% after spending quite some time with Phy checking on the results and understanding the merge and split functions, 
% come back here and run Kilosort's automated merging strategy. This block
% will overwrite the previous results and python files. Load the results in
% Phy again: there should be no merges left to do (with the default simulation), but perhaps a few splits
% / cleanup. On realistic data (i.e. not this simulation) there will be drift also, which will usually
% mean there are merges left to do even after this step. 
% Kilosort's AUTO merges should not be confused with the "best" merges done inside the
% benchmark (those are using the real ground truth!!!)

% LFP
set(0,'DefaultFigureWindowStyle','normal')
LFP = fastpreprocess_filtering(flip(Intan.allIntan,1),8192); %Only run for PFF data
% LFP = fastpreprocess_filtering(Intan.allIntan,8192);

LFP = bestLFP(LFP);
LFP = bandFilter(LFP,'depth'); % Extract LFPs based on 'depth' or 'single'
% LFPplot(LFP)
% %% Beta Band Analysis
% LFP = betaBurstDetection(LFP);
% %% CSD
% [CSDoutput]  = CSD(LFPavg'/1E6,8192,2E-5);
% Looking at single units
% set(0,'DefaultFigureWindowStyle','docked')
Spikes = singleUnitAnalysis(fpath,VR_data); % VR_data.Time{1} = data(:,2); VR_data.Position{1} = data(:,1);
% Calculate Depth profile
set(0,'DefaultFigureWindowStyle','normal')
load chanMap % use for PFF
% load UCLA_chanMap_fixed
[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] =...
    spikeTemplatePosition(fpath,ycoords);
figure,
plot(waveforms(:,20:end-20)','color',[0.5 0.5 0.5 0.25]), hold on;
% figure,plot(mean(waveforms(1:30,20:end-20),1),'k','LineWidth',2)
% test

Spikes = spikeDepthPlot(Spikes,templateDepths);

%%
% Time-Frequency Analysis
[TimeFreq,LFP,betaGroup,Spikes] = tfAnalysis(Spikes,LFP,1); %Behavior state running 1 (0 rest)
[TimeFreq,LFP,betaGroupRest,Spikes] = tfAnalysis(Spikes,LFP,0,TimeFreq); %Behavior state running 1 (0 rest)

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
itpcstats = tfStats(TimeFreq);ylim([0 0.4])

savepath = 'Y:\Hammad\Ephys\PFFProject\RebuttelITPC\';
sessionName = [savepath,filename,'_itpcStats.mat'];
% save(sessionName,"IntanBehaviour","fpath","parameters","-v7.3");
save(sessionName,"pathname","itpcstats","Spikes","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
disp('Data Saved!')
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
%  pxSpecs = statePowerSpec(TimeFreq,LFP);  old version
idx = find(s.sorted_probe_wiring(:,2)==0);
params.Fs = 1024;
params.fpass = [1 100];
params.tapers = [5 9];
movingwin = [0.5 0.5];
params.pad = 1;
params.trialave = 1;
params.err = [2 0.05];
[Srun,f,SerrRun] = mtspectrumc(squeeze(TimeFreq.tfRun.fullLFP(:,idx(end),:)),params);
[Srest,f,SerrRest] = mtspectrumc(squeeze(TimeFreq.tfRest.fullLFP(:,idx(end),1:50)),params);
figure,semilogy(f,Srun);hold on;semilogy(f,SerrRun);title('Run'),ylim([0 20E2]),box off,set(gca,'TickDir','out'),set(gca,'FontSize',14)
figure,semilogy(f,Srest);hold on;semilogy(f,SerrRest);title('Rest'),ylim([0 20E2]),box off,set(gca,'TickDir','out'),set(gca,'FontSize',14)
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
[L23RunAvg,L5RunAvg] = spikeRaster(Spikes,1,5); %Spikes, flag (0/1) for rest/run, and scatter size (sz)
[L23RestAvg,L5RestAvg] = spikeRaster(Spikes,0,5); %Spikes, flag (0/1) for rest/run
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
end
for i = 1:size(betaGroup,2) % Checks electrode size for median
    disp(['Electrode: ' num2str(i)])
    [peakAlign{i},csd{i},betaNorm{i},f,bstats(i)] = betaAnalysis(betaGroup(i).electrode);
end 

% Beta Event Rate
betaEventRate = [];
betaEventDuration = [];
for i = 1:64
    betaEventRate = [betaEventRate; bstats(i).betaER];
    betaEventDuration= [betaEventDuration;bstats(i).betaDuration];
end
figure,boxplot(betaEventRate,'Plotstyle','compact'),ylim([1 6]),box off
figure,boxplot(betaEventDuration,'Plotstyle','compact'),box off
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
betaStats(bstats,LFPdepth)

%% Plot beta traces for each electrode
figure, hold on
for i = 1:size(peakAlign,2) % Checks electrode size for median
    mPeakAlign(:,i) = mean(peakAlign{i},1);
    plot(peakAlign{i}')
end
mPeakAlign = mPeakAlign(:,[1:43,45:64]);
figure,stack_plot(flip(mPeakAlign'),0.2,1.5)
normPeakAlign = (mPeakAlign-min(mPeakAlign,[],'all'))/(max(mPeakAlign,[],'all')-min(mPeakAlign,[],'all'));
figure,imagesc(0:250,LFPdepth,interp2(smoothdata(normPeakAlign'))),caxis([0.1 1])
%% Plot beta CSD for each electrode
for i = 1:size(csd,2) % Checks electrode size for median
    mcsd(:,:,i) = mean(csd{i},3);
end
mcsd = mean(mcsd,3);
figure,imagesc(0:250,LFPdepth,interp2(smoothdata((mcsd')),2)),colormap(jet),caxis([-0.2 0.2])

%%
behaviorflag = 1;
spikeTriggeredBeta = betaEventPSH(betaGroup,Spikes,behaviorflag); %set behavior flag 0 or 1 for rest/run
%% Layer specific spike rate
% Convert cell x trial to single array (cos we don't care what trial this
% all happens on
L23RunSR = Spikes.spikeRate.L23RunSR;
L23RunSR = L23RunSR(:); %Reshape to 1d
L23RunSR(L23RunSR<1) = []; %Remove zeros

L4RunSR = Spikes.spikeRate.L4RunSR;
L4RunSR = L4RunSR(:); %Reshape to 1d
L4RunSR(L4RunSR<1) = []; %Remove zeros

L5RunSR = Spikes.spikeRate.L5RunSR;
L5RunSR = L5RunSR(:); %Reshape to 1d
L5RunSR(L5RunSR<1) = []; %Remove zeros

L23RestSR = Spikes.spikeRate.L23RestSR;
L23RestSR = L23RestSR(:); %Reshape to 1d
L23RestSR(L23RestSR<1) = [];

L4RestSR = Spikes.spikeRate.L4RestSR;
L4RestSR = L4RestSR(:); %Reshape to 1d
L4RestSR(L4RestSR<1) = [];

L5RestSR = Spikes.spikeRate.L5RestSR;
L5RestSR = L5RestSR(:); %Reshape to 1d
L5RestSR(L5RestSR<1) = [];

%Avg SR
avgSR = vertcat(Spikes.spikeRate.L23avgSR,Spikes.spikeRate.L4avgSR,Spikes.spikeRate.L5avgSR);
xAvgSR = [repmat({'L23'},length(Spikes.spikeRate.L23avgSR),1);repmat({'L4'},length(Spikes.spikeRate.L4avgSR),1);repmat({'L5'},length(Spikes.spikeRate.L5avgSR),1)];
figure,boxplot(avgSR,xAvgSR,'plotstyle','compact'),title('Average SR'), hold on
set(gca, 'YScale', 'log')
ylim([00.1 100])
x = [1.1*ones(length(Spikes.spikeRate.L23avgSR),1);2.1*ones(length(Spikes.spikeRate.L4avgSR),1);3.1*ones(length(Spikes.spikeRate.L5avgSR),1)];
scatter(x,avgSR,'filled','k')

%Rest SR
restSR = vertcat(L23RestSR,L4RestSR,L5RestSR);
xRestSR = [repmat({'L23'},length(L23RestSR),1);repmat({'L4'},length(L4RestSR),1);repmat({'L5'},length(L5RestSR),1)];
figure,boxplot(restSR,xRestSR,'plotstyle','compact'),title('Resting SR'),hold on
ax = gca;
ax.YAxis.Scale ='log';
ylim([00.1 100])
x = [1.1*ones(length(L23RestSR),1);2.1*ones(length(L4RestSR),1);3.1*ones(length(L5RestSR),1)];
scatter(x,restSR,'filled','k')

% Run SR
runSR = vertcat(L23RunSR,L4RunSR,L5RunSR);
xRunSR = [repmat({'L23'},length(L23RunSR),1);repmat({'L4'},length(L4RunSR),1);repmat({'L5'},length(L5RunSR),1)];
figure,boxplot(runSR,xRunSR,'plotstyle','compact'),title('Running SR'),hold on
ax = gca;
ax.YAxis.Scale ='log';
ylim([00.1 100])
x = [1.1*ones(length(L23RunSR),1);2.1*ones(length(L4RunSR),1);3.1*ones(length(L5RunSR),1)];
scatter(x,runSR,'filled','k')
%%
%Histogram for each layer by brain state
[~,edges] = histcounts(log10(L23RunSR));
figure, histogram(L23RunSR,10.^edges),title('Layer 2/3'), hold on
histogram(L23RestSR(1:length(L23RunSR),1),25.^edges);
set(gca, 'xscale','log')
xlim([00.1 100])                 

[~,edges] = histcounts(log10(L4RunSR));
figure, histogram(L4RunSR,10.^edges),title('Layer 4'), hold on
histogram(L4RestSR(1:length(L4RunSR),1),25.^edges);
set(gca, 'xscale','log')
xlim([00.1 100])

[~,edges] = histcounts(log10(L5RunSR));
figure, histogram(L5RunSR,10.^edges),title('Layer 5'), hold on
histogram(L5RestSR(1:length(L5RunSR),1),25.^edges);
set(gca, 'xscale','log')
xlim([00.1 100])

%% histogram stats for beta event
%
figure,histogram(bstats(1:64,1),50:100),title('Event Duration')
figure,histogram(bstats(1:64,3),1:0.1:5),title('Event Rate')
figure,
%% Beta states
% Number of beta per electrode
for i = 1:size(betaGroup,2)
    runBetaNum(i) = mean(betaGroup(i).electrode.betaBurst.NumDetectedBeta);
end
% Number of beta per trial per electrode
for i = 1:size(betaGroup,2)
    runBetaDepth(i,:) = betaGroup(i).electrode.betaBurst.NumDetectedBeta;
end

%Number of beta per electrode during rest
for i = 1:size(betaGroupRest,2)
    restBetaNum(i) = mean(betaGroupRest(i).electrode.betaBurst.NumDetectedBeta);
end

% Number of beta per trial per electrode
for i = 1:size(betaGroupRest,2)
    restBetaDepth(i,:) = betaGroupRest(i).electrode.betaBurst.NumDetectedBeta;
end
%% 
figure,
colorH = jet(4);
for i = 1:64
    if i<=30
        col = colorH(1,:);
    else
        col = colorH(4,:);
    end
    plot(bstats(i).freqPks,'.','color',col), hold on
    freqPks = bstats(i).freqPks;
    lowBeta{i} = freqPks(freqPks<=21);
    highBeta{i} = freqPks(freqPks>21);
end
L23lowBeta = lowBeta(1:30);
L23lowBeta = horzcat(L23lowBeta{:});
L5lowBeta = lowBeta(31:64);
L5lowBeta = horzcat(L5lowBeta{:});

L23highBeta = highBeta(1:30);
L23highBeta = horzcat(L23highBeta{:});
L5highBeta = highBeta(31:64);
L5highBeta = horzcat(L5highBeta{:});

figure,histogram(L23highBeta,20:1:30),title('L23 High Beta')
figure,histogram(L23lowBeta,13:1:20),title('L23 Low Beta')
figure,histogram(L5highBeta,20:1:30),title('L5 High Beta')
figure,histogram(L5lowBeta,13:1:20),title('L5 Low Beta')


%%
figure,histogram(restBetaNum,1:0.25:10),hold on
histogram(runBetaNum,1:0.25:10)
figure,lineError(1:64,restBetaDepth')
figure,lineError(1:64,runBetaDepth')
%%
bothBeta = [restBetaNum;runBetaNum];
figure,bar(mean(bothBeta,2))
hold on
er = errorbar(mean(bothBeta,2),std(bothBeta'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
