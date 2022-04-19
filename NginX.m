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
ParpoolConfig
fpath    = Intan.path; % where on disk do you want the analysis? ideally and SSD...
pathToYourConfigFile = strcat(pwd,'/main/'); % for this example it's ok to leave this path inside the repo, but for your own config file you *must* put it somewhere else!  
run(fullfile(pathToYourConfigFile, 'config_eMouse.m'))
make_UCLAMouseChannelMap(fpath); % Creates channel map for electrode array
%%
% filtData = preprocess_filtering(Intan.allIntan(:,1:400000),Intan.t_amplifier);
%% Ripples
% [Ripples,filtData] = SWR(Intan);
% plotClusterless(Ripples,filtData,Intan)


%% Kilosort Analysis
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

% LFP
set(0,'DefaultFigureWindowStyle','normal')
LFP = fastpreprocess_filtering(flip(Intan.allIntan,1),8192);
LFP = bestLFP(LFP);
LFP = bandFilter(LFP,'depth'); % Extract LFPs based on 'depth' or 'single'
% LFPplot(LFP)
% %% Beta Band Analysis
% LFP = betaBurstDetection(LFP);
% %% CSD
% [CSDoutput]  = CSD(flip(LFP.LFP(:,1:1024)'/1E6,1),1024,2E-5);
% Looking at single units
% set(0,'DefaultFigureWindowStyle','docked')
Spikes = singleUnitAnalysis(fpath,VR_data);
% Calculate Depth profile
set(0,'DefaultFigureWindowStyle','normal')
load chanMap
[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] =...
    spikeTemplatePosition(fpath,ycoords);
% figure,
% plot(waveforms(:,20:end-20)','color',[0.5 0.5 0.5 0.25]), hold on;
% figure,plot(mean(waveforms(1:30,20:end-20),1),'k','LineWidth',2)
% test

Spikes = spikeDepthPlot(Spikes,templateDepths);
% Time-Frequency Analysis
[TimeFreq,LFP,betaGroup,Spikes] = tfAnalysis(Spikes,LFP,1); %Behavior state running 1 (0 rest)
[TimeFreq,LFP,betaGroupRest,Spikes] = tfAnalysis(Spikes,LFP,0); %Behavior state running 1 (0 rest)
%%
% [TimeFreq,LFP,betaGroupRest,Spikes] = tfAnalysis(Spikes,LFP,0,TimeFreq); %Behavior state running 1 (0 rest)

% plotTF(TimeFreq,LFP)
% TF stats of depth
% TimeFreq.tf = TimeFreq.tfRun;
% stats = tfStats(TimeFreq);
% tfDepth = TimeFreq.tf.depth;
betaGammaCoupling = gammaBetaCoupling(TimeFreq.gamma,TimeFreq.beta);
betaGammam = mean(betaGammaCoupling,3);
figure,imagesc(-179:20:180,1:64,interp2(betaGammam')),colormap(jet)
figure,plot(mean(betaGammam))
% figure,imagesc(betaGammam(:,1:40)'),colormap(jet)
% figure,imagesc(betaGammam(:,41:64)'),colormap(jet)

%% 
spikeRaster(Spikes)
%% Beta Analysis for each electrode
for i = 1:size(LFP.medianLFP,1) % Checks electrode size for median
    disp(['Electrode: ' num2str(i)])
    [peakAlign{i},csd{i},betaNorm{i},f,bstats(i)] = betaAnalysis(betaGroup(i).electrode,LFP.LFP);
end 
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
    bstats = cell2mat(struct2cell(bstats));
    bstats = squeeze(bstats)';
end
figure('Name','Beta Duration'),bar(LFPdepth,bstats(:,1),'BarWidth',1),set(gcf, 'Position',  [100, 100, 500, 500])
figure('Name','Beta Amplitude'),bar(LFPdepth,bstats(:,2),'BarWidth',1),set(gcf, 'Position',  [100, 100, 500, 500])
figure('Name','Beta Event Rate'),bar(LFPdepth,bstats(:,3),'BarWidth',1),set(gcf, 'Position',  [100, 100, 500, 500])
% Bar plot for each layer
betaStats(bstats,LFPdepth)

%% Plot beta traces for each electrode
figure, hold on
for i = 1:size(peakAlign,2) % Checks electrode size for median
    mPeakAlign(:,i) = mean(peakAlign{i},1);
    plot(peakAlign{i}')
end
figure,stack_plot(flip(mPeakAlign'),0.2,1.5)
normPeakAlign = (mPeakAlign-min(mPeakAlign,[],'all'))/(max(mPeakAlign,[],'all')-min(mPeakAlign,[],'all'));
figure,imagesc(0:250,LFPdepth,smoothdata(normPeakAlign')),caxis([0 1])
%% Plot beta CSD for each electrode
for i = 1:size(csd,2) % Checks electrode size for median
    mcsd(:,:,i) = mean(csd{i},3);
end
mcsd = mean(mcsd,3);
figure,imagesc(0:250,LFPdepth,interp2(smoothdata((mcsd')),2)),colormap(jet),caxis([-0.2 0.2])

%% Beta Event Triggered Spike Rate
% temp window
for i = 1:11 % trials
    for ii = 1:30 % L2/3 beta events
        betaEventTemp = cell2mat(betaGroup(ii).electrode.betaBurst.detectedBeta(i));
        for j = 1:size(betaEventTemp,1)
            for jj = 1:size(Spikes.spikes.L23Run(i).Cell,2)
                spikeInBeta{jj,j,ii,i} = find(Spikes.spikes.L23Run(i).Cell{jj}>=(betaEventTemp(j,2)-.075)...
                    & Spikes.spikes.L23Run(i).Cell{jj}<=((betaEventTemp(j,2)+.075)));
            end
        end
    end
end
%% L2/3
for i = 1:11 % trials
    count = 1;
    for ii = 1:30 % L2/3 beta events
        betaEventTemp = cell2mat(betaGroup(ii).electrode.betaBurst.detectedBeta(i));
        for j = 1:size(betaEventTemp,1)
            spikeTemp = vertcat(Spikes.spikes.L23Run(i).Cell{:});
            spikeInBetaTemp = find(spikeTemp>=(betaEventTemp(j,1)-betaGroup(ii).electrode.betaBurst.window(i,1))...
                & spikeTemp<=(betaEventTemp(j,3)-betaGroup(ii).electrode.betaBurst.window(i,1)));
            spikeTriggeredBeta(i).L23.spike{count} = spikeTemp(spikeInBetaTemp);
            spikeTriggeredBeta(i).L23.betaEvent(count,:) = betaEventTemp(j,:);
            
            start = spikeTriggeredBeta(i).L23.betaEvent(count,1)*1024;
            stop = spikeTriggeredBeta(i).L23.betaEvent(count,3)*1024;
            betaLine = betaGroup(ii).electrode.beta_band(start:stop);
%             plot(spikeTemp(spikeInBetaTemp),i*ones(length(spikeInBetaTemp)),'.'),hold on
            plot(betaLine),hold on
            spikeTriggeredBeta(i).L23.betaLFP{count} = betaLine;
            count = count+1;
        end
    end
end
%% Layer 5
for i = 1:11 % trials
    count = 1;
    for ii = 30:64 % L5beta events
        betaEventTemp = cell2mat(betaGroup(ii).electrode.betaBurst.detectedBeta(i));
        for j = 1:size(betaEventTemp,1)
            spikeTemp = vertcat(Spikes.spikes.L23Run(i).Cell{:});
            spikeInBetaTemp = find(spikeTemp>=(betaEventTemp(j,1)-betaGroup(ii).electrode.betaBurst.window(i,1))...
                & spikeTemp<=(betaEventTemp(j,3)-betaGroup(ii).electrode.betaBurst.window(i,1)));
            spikeTriggeredBeta(i).L5.spike{count} = spikeTemp(spikeInBetaTemp);
            spikeTriggeredBeta(i).L5.betaEvent(count,:) = betaEventTemp(j,:);
            
            start = spikeTriggeredBeta(i).L5.betaEvent(count,1)*1024;
            stop = spikeTriggeredBeta(i).L5.betaEvent(count,3)*1024;
            betaLine = betaGroup(ii).electrode.beta_band(start:stop);
            plot(betaLine),hold on
            spikeTriggeredBeta(i).L5.betaLFP{count} = betaLine;
            count = count+1;
        end
    end
end
%%
betaEventPSH(betaGroup,Spikes)
%%
figure,
for i = 1:35
subplot(2,1,1),plot(spikeTriggeredBeta(1).L23.betaLFP{i}),hold on
subplot(2,1,2),plot(spikeTriggeredBeta(1).L23.spike{i},i*ones(length(spikeTriggeredBeta(1).L23.spike{i})),'.'),hold on

end

%%
for i = 1:35

end
%%
    for j = 30:64 %L5 beta events
        betaGroup(j).electrode.betaBurst.detectedBeta(i);
    end
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

%% histogram stats
figure,histogram(bstats(1:64,1),50:100),title('Event Duration')
figure,histogram(bstats(1:64,3),1:0.25:10),title('Event Rate')
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
