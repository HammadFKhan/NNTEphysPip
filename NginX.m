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
LFP = fastpreprocess_filtering(Intan.allIntan,8192);
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
%%
[TimeFreq,LFP,betaGroupRest,Spikes] = tfAnalysis(Spikes,LFP,0); %Behavior state running 1 (0 rest)

% [TimeFreq,LFP,betaGroupRest,Spikes] = tfAnalysis(Spikes,LFP,0,TimeFreq); %Behavior state running 1 (0 rest)

% plotTF(TimeFreq,LFP)
% TF stats of depth
% TimeFreq.tf = TimeFreq.tfRun;
% stats = tfStats(TimeFreq);
% tfDepth = TimeFreq.tf.depth;
betaGammaCoupling = gammaBetaCoupling(TimeFreq.gamma,TimeFreq.beta);
betaGammam = mean(betaGammaCoupling,3);
figure,imagesc(betaGammam'),colormap(jet)
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
for i = 1:size(peakAlign,2) % Checks electrode size for median
    mPeakAlign(:,i) = mean(peakAlign{i},1);
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

%% Construct ERD/ERS metric
figure,hold on
L23EventISI = [];
L23RunEventISI = [];
L23InitiateEventISI = [];
%L2/3
for ii = 30:64
    for i = 1:size(betaGroup(ii).electrode.betaBurst.NumDetectedBeta,1)
        if isempty(betaGroup(ii).electrode.betaBurst.detectedBeta{i})
            t = 0;t1 = 0;
        else
            t = betaGroup(ii).electrode.betaBurst.detectedBeta{i}(:,2);
            t1 = betaGroup(ii).electrode.betaBurst.window(i,1);
            win = t-t1;
            plot(win,i*ones(1,length(t)),'.')
            % ERD (Event Rate Desynchonization)
            %         A = win(win<=1);
            %         B = win(win>1);
            %         ERD(i) = (length(A)-length(B))/(length(A)+length(B));
            %Beta Event Interval
            L23EventISI{i} = abs(diff(win));
            L23RunEventISI{i} = abs(diff(win(win>=1))); % Run Window
            L23InitiateEventISI{i} = abs(diff(win(win<1))); % Initiate Window
            if isempty(L23RunEventISI{i}) && isempty(L23InitiateEventISI{i})
                L23RunEventISI{i} = -.5;
                L23InitiateEventISI{i} = -.5;
            elseif isempty(L23RunEventISI{i}) 
                L23RunEventISI{i} = -.5;
            elseif isempty(L23InitiateEventISI{i})
                L23InitiateEventISI{i} = -.5;
            end
        end
    end
end

L23EventISI = vertcat(L23EventISI{:});
L23InitiateEventISI = vertcat(L23InitiateEventISI{:});
L23RunEventISI = vertcat(L23RunEventISI{:});

figure('Name','Layer 2/3'),histogram(L23EventISI,0:.1:2)
figure('Name','Initiate Layer 2/3'),yyaxis left,histogram(L23InitiateEventISI,0:.1:2),title('Initiate Layer2/3'),ylim([0 45])
yyaxis right,histogram(L23InitiateEventISI(L23InitiateEventISI==-0.5),-1:0.1:0),ylim([0 350])
figure('Name','Run Layer 2/3'),yyaxis left,histogram(L23RunEventISI,0:.1:2),title('Run Layer 2/3')
yyaxis right,histogram(L23RunEventISI(L23RunEventISI==-0.5),-1:0.1:0),ylim([0 350])


% L5
figure,hold on
L5EventISI = [];
L5InitiateEventISI = [];
L5RunEventISI = [];
for ii = 1:30
    for i = 1:size(betaGroup(ii).electrode.betaBurst.NumDetectedBeta,1)
        if isempty(betaGroup(ii).electrode.betaBurst.detectedBeta{i})
            t = 0;t1 = 0;
        else
            t = betaGroup(ii).electrode.betaBurst.detectedBeta{i}(:,2);
            t1 = betaGroup(ii).electrode.betaBurst.window(i,1);
            win = t-t1;
            plot(win,i*ones(1,length(t)),'.')
            % ERD (Event Rate Desynchonization)
            %         A = win(win<=1);
            %         B = win(win>1);
            %         ERD(i) = (length(A)-length(B))/(length(A)+length(B));
            L5EventISI{i} = abs(diff(win));
            L5RunEventISI{i} = abs(diff(win(win>=1))); % Run Window
            L5InitiateEventISI{i} = abs(diff(win(win<1))); % Initiate Window
            if isempty(L5RunEventISI{i}) && isempty(L5InitiateEventISI{i})
                L5RunEventISI{i} = -.5;
                L5InitiateEventISI{i} = -.5;
            elseif isempty(L5RunEventISI{i})
                L5RunEventISI{i} = -.5;
            elseif isempty(L5InitiateEventISI{i})
                L5InitiateEventISI{i} = -.5;
            end
        end
    end
end
L5EventISI = vertcat(L5EventISI{:});
L5InitiateEventISI = vertcat(L5InitiateEventISI{:});
L5RunEventISI = vertcat(L5RunEventISI{:});

figure('Name','Layer 5'),histogram(L5EventISI,0:.1:2)
figure('Name','Initiate Layer 5'),yyaxis left,histogram(L5InitiateEventISI,0:.1:2),title('Initiate Layer 5'),ylim([0 45])
yyaxis right,histogram(L5InitiateEventISI(L5InitiateEventISI==-0.5),-1:0.1:0),ylim([0 350])
figure('Name','Run Layer 5'),yyaxis left,histogram(L5RunEventISI,0:.1:2),title('Run Layer 5')
yyaxis right,histogram(L5RunEventISI(L5RunEventISI==-0.5),-1:0.1:0),ylim([0 350])
% Rest State
figure,hold on
L23RestEventISI = [];
%L2/3
for ii = 30:64
    for i = 1:size(betaGroupRest(ii).electrode.betaBurst.NumDetectedBeta,1) % Keep the same number of trials
        if isempty(betaGroupRest(ii).electrode.betaBurst.detectedBeta{i})
            t = 0;t1 = 0;
        else
            t = betaGroupRest(ii).electrode.betaBurst.detectedBeta{i}(:,2);
            t1 = betaGroupRest(ii).electrode.betaBurst.window(i,1);
            win = t-t1;
            win = win(win<1); %% split in half to make the inwdow the same for run and intitate
            plot(win,i*ones(1,length(win)),'.')
            %Beta Event Interval
            L23RestEventISI{i} = abs(diff(win));
            if isempty(L23RestEventISI{i})
                L23RestEventISI{i} = -.5;
            end
        end
    end
end
L23RestEventISI = vertcat(L23RestEventISI{:});
figure('Name','Rest Layer 2/3'),yyaxis left,histogram(L23RestEventISI,0:.1:2),title('Rest Layer 2/3')
yyaxis right,histogram(L23RestEventISI(L23RestEventISI==-0.5),-1:0.1:0),ylim([0 350])
% L5
figure,hold on
L5RestEventISI = [];
for ii = 1:30
    for i = 1:size(betaGroupRest(ii).electrode.betaBurst.NumDetectedBeta,1)
        if isempty(betaGroupRest(ii).electrode.betaBurst.detectedBeta{i})
            t = 0;t1 = 0;
        else
            t = betaGroupRest(ii).electrode.betaBurst.detectedBeta{i}(:,2);
            t1 = betaGroupRest(ii).electrode.betaBurst.window(i,1);
            win = t-t1;
            win = win(win<1); %% split in half to make the inwdow the same for run and intitate
            plot(win,i*ones(1,length(win)),'.')
            L5RestEventISI{i} = abs(diff(win));
            if isempty(L5RestEventISI{i})
                L5RestEventISI{i} = -.5;
            end
        end
    end
end
L5RestEventISI = vertcat(L5RestEventISI{:});
figure('Name','Rest Layer 5'),yyaxis left,histogram(L5RestEventISI,0:.1:2),title('Rest Layer 5')
yyaxis right,histogram(L5RestEventISI(L5RestEventISI==-0.5),-1:0.1:0),ylim([0 350])
%%
figure,boxplot(L23InitiateEventISI(L23InitiateEventISI>0),'PlotStyle','compact'),ylim([0 1]),title('L23 Initiate')
figure,boxplot(L23RunEventISI(L23RunEventISI>0),'PlotStyle','compact'),ylim([0 1]),title('L23 Run')
figure,boxplot(L5InitiateEventISI(L5InitiateEventISI>0),'PlotStyle','compact'),ylim([0 1]),title('L5 Initiate')
figure,boxplot(L5RunEventISI(L5RunEventISI>0),'PlotStyle','compact'),ylim([0 1]),title('L5 Run')
%%
l23run = [length(L23RunEventISI(L23RunEventISI>0));length(L23RunEventISI(L23RunEventISI==-0.5))]/size(L23RunEventISI,1);
figure,pie(l23run),title('L23 Run')
legend('Burst','Single')

l23init = ([length(L23InitiateEventISI(L23InitiateEventISI>0));length(L23InitiateEventISI(L23InitiateEventISI==-0.5))])/size(L23InitiateEventISI,1);
figure,pie(l23init),title('L23 Initiate')
legend('Burst','Single')

l5run = [length(L5RunEventISI(L5RunEventISI>0));length(L5RunEventISI(L5RunEventISI==-0.5))]/size(L5RunEventISI,1);
figure,pie(l5run),title('L5 Run')
legend('Burst','Single')

l5init = [length(L5InitiateEventISI(L5InitiateEventISI>0));length(L5InitiateEventISI(L5InitiateEventISI==-0.5))]/size(L5InitiateEventISI,1);
figure,pie(l5init),title('L5 Initiate')
legend('Burst','Single')

l23rest = [length(L23RestEventISI(L23RestEventISI>0));length(L23RestEventISI(L23RestEventISI==-0.5))]/size(L23RestEventISI,1);
figure,pie(l23rest),title('L2/3 Rest')
legend('Burst','Single')

l5rest = [length(L5RestEventISI(L5RestEventISI>0));length(L5RestEventISI(L5RestEventISI==-0.5))]/size(L5RestEventISI,1);
figure,pie(l5rest),title('L5 Rest')
legend('Burst','Single')
%%
%%
figure,boxplot(L23EventISI,'PlotStyle','compact'),ylim([0 2]),title('L23')
figure,boxplot(L5EventISI,'PlotStyle','compact'),ylim([0 2]),title('L5')
figure,boxplot(L23RestEventISI,'PlotStyle','compact'),ylim([0 2]),title('L23 Rest')
figure,boxplot(L5RestEventISI,'PlotStyle','compact'),ylim([0 2]),title('L5 Rest')
%%
figure,bar(sort(ERD)),ylim([-1 1])
figure,bar(sort(ERDrest)),ylim([-1 1])
figure, histogram(ERDrest,-1:0.5:1),hold on
histogram(ERD,-1:0.5:1)
% ERD probability
ERDprop = length(ERD(abs(ERD)>=0.6))/length(ERD);
ERDpropRest = length(ERDrest(abs(ERDrest)>=0.6))/length(ERDrest);
figure,bar([ERDpropRest;ERDprop]),hold on
% er = errorbar([ERDpropRest;ERDprop],std([ERDpropRest;ERDprop]'));    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  

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
