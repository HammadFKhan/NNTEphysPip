load Y:\Austin\ShankData\240612_shank_eOPN_M1\loadme
load Y:\Austin\ShankData\240612_shank_eOPN_M1\IntanBehaviour

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
%load chanMap64F2
load chanMap64Sharp
[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms, max_site] =...
    spikeTemplatePosition(data.fpath,ycoords,[]); % 'invert'
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
%% Calculate trial PSTH for passive stimulation via pole

% TODO: Make a parameters structure to track experimental settings. Need
% this for correct analysis
parameters.windowBeforePole = 0.5;
parameters.windowAfterPole = 0.75;
parameters.experiment = 'passive'; % passive - pole, active - active touch
parameters.opto = 0; % 1 - opto ON , 0 - opto OFF
parameters.Fs = 1000; % Eventual downsampled data
parameters.ts = 1/parameters.Fs;

Spikes = polePSTH(Spikes,IntanBehaviour);
%%% save spike output data to load into gui
savepath = fullfile(path,['spks4sorting','.mat']);
save(savepath,'Spikes','-v7.3')
%% Basic spike analysis
% z-score spike rates
if exist('parameters','var')
    IntanBehaviour.parameters = parameters;
end
if exist('goodSpkComponents','var')
    Spikes.goodSpkComponents = unique(goodSpkComponents);
else 
    Spikes.goodSpkComponents = 1:length(Spikes.Clusters);
end
Spikes = rejectSpikespassive(Spikes,0.60,0.25,parameters); % Reject spikes here for further analysis: (Spikes,fractionTrials,cutoffFR,parameters)
[fpath,name,exts] = fileparts(ds_filename);
sessionName = [fpath,'/','Spikes.mat'];
% save(sessionName,"IntanBehaviour","fpath","parameters","-v7.3");
save(sessionName,"Spikes","IntanBehaviour","fpath","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
disp('Saved!')

%% Plot some rep units
figure,hold on
time = -parameters.windowBeforePole*parameters.Fs:parameters.windowAfterPole*parameters.Fs;
for n = 54
subplot(2,1,[1]),Show_Spikes(Spikes.PSTH.pole.spks{n}),axis off
subplot(2,1,[2]),bar(time,smoothdata(Spikes.PSTH.pole.spkRates(n,:)),'FaceColor',[28/255 117/255 188/255],'EdgeColor','none')
axis tight, box off, set(gca,'TickDir','out')
set(gca,'fontsize',16)
xline(0)
end
%% plot heatmap of spikes arranged by depth
figure,imagesc(time,1:size(Spikes.PSTH.pole.spkRates,1),Spikes.PSTH.pole.spkRates);colormap(hot),caxis([0 50]),colorbar
figure,plot(mean(Spikes.PSTH.pole.spkRates)),box off,hold on
plot(mean(Spikes.PSTH.pole.spkRates)+std(Spikes.PSTH.pole.spkRates)/sqrt(size(Spikes.PSTH.pole.spkRates,1)));
plot(mean(Spikes.PSTH.pole.spkRates)-std(Spikes.PSTH.pole.spkRates)/sqrt(size(Spikes.PSTH.pole.spkRates,1)))
%% eOPN analysis
baseline = 1:112;
eOPN = 113:212;
baselineFR = cellfun(@(x) mean(x(baseline,:)),Spikes.PSTH.pole.spks,'UniformOutput',false);
baselineFR = smoothdata(vertcat(baselineFR{:})*1000,2,'gaussian',20);
eOPNFR = cellfun(@(x) mean(x(eOPN,:)),Spikes.PSTH.pole.spks,'UniformOutput',false);
eOPNFR = smoothdata(vertcat(eOPNFR{:})*1000,2,'gaussian',20);
figure,plot(time,mean(baselineFR)),box off, hold on
plot(time,mean(baselineFR)+std(baselineFR)/sqrt(size(baselineFR,1)));
plot(time,mean(baselineFR)-std(baselineFR)/sqrt(size(baselineFR,1)));

figure,plot(time,mean(eOPNFR)),box off, hold on
plot(time,mean(eOPNFR)+std(eOPNFR)/sqrt(size(eOPNFR,1)));
plot(time,mean(eOPNFR)-std(eOPNFR)/sqrt(size(eOPNFR,1)));
%% Prep data for VarMean analysis
S1Data = struct();
count = 1;

for n = 1:length(Spikes.PSTH.pole.spks)
    S1Data(count).spikes = logical(Spikes.PSTH.pole.spks{n}(baseline,:));
    count = count+1;
end

%% Calculate Neural variability via mean-matched FF
% See: Churchland MM et al. (2010)  Stimulus onset quenches neural variability: a widespread cortical phenomenon. Nat. Neurosci.

% DEFINE Data type here
% Call FF function and calculate
M1Data = S1Data;
ResultsBaseline = calcNeuralVariance(M1Data);
M1Data = S1DataeOPN;
ResultseOPN = calcNeuralVariance(M1Data);

%% Check Fakerized data 
% This is critical!! Here we fakerize our actual data to check if the effect
% is real. 
% The result is 'fake' data whose basic properties (e.g. mean firing rate) are identical to the
% original data, but where every trial has the same underlying rate, and the only across-trial
% variability is due to spiking statistics. 

% The function works by redistributing the spikes in a given ms across all trials for that
% neuron/condition.  Thus, any within-trial spiking autocorrelation is removed.  The result
% is poisson spiking statistics (within 1 ms resolution) with no change in the mean rate.

%
% Using default Poisson statistics, the Fano factor for 'Fakerized' data should be very close to 1. 
% It will be slightly less due to the 1 ms refractory period imposed by the data format, and will 
% drop slightly if the firing rate rises.  However, when matching spike-count distributions 
% (using VarVsMean the latter effect should dissapear.  If the Fano Factor produced by VarVsMean does 
% NOT remain constant (probably near 0.95)when using Fakerized data (and when matching dists) then 
% there must be an artifact in the analysis. 
% I statistically plot this out in the next part :)

M1Data = S1Data;
FakResulthit = calcFakNeuralVariace(M1Data);
M1Data = M1Datamiss;
FakResultmiss = calcFakNeuralVariace(M1Data);
M1Data = M1DataMIFA;
FakResultFA = calcFakNeuralVariace(M1Data);
%% Calc FF stats and make structure for further analysis
Stats.Baseline = calcStats(ResultsBaseline);
Stats.eOPN = calcStats(ResultseOPN);
%%

Stats.FakHit = calcStats(FakResulthit);
Stats.FakMiss = calcStats(FakResultmiss);
Stats.FakFA = calcStats(FakResultFA);
%% Plot it out
dat1 = Stats.Baseline.FFdrop;
dat2 = Stats.eOPN.FFdrop;

% dat1 = Stats.FakHit.FFdrop;
% dat2 = Stats.FakMiss.FFdrop;
% dat3 = Stats.FakFA.FFdrop;

temp = zeros(max([length(dat1) length(dat2)]),2);
temp(1:length(dat1),1) = dat1;
temp(1:length(dat2),2) = dat2;
figure
subplot(121),customBoxplot(temp),set(gca,'tickdir','out','fontsize',16),box off,ylim([-1 1.5])

% dat1 = Stats.FakBaseline.FFstim;
% dat2 = Stats.FakeOPN.FFstim;

dat1 = Stats.Baseline.FFstim;
dat2 = Stats.eOPN.FFstim;

temp = zeros(max([length(dat1) length(dat2)]),2);
temp(1:length(dat1),1) = dat1;
temp(1:length(dat2),2) = dat2;
subplot(122),customBoxplot(temp),set(gca,'tickdir','out','fontsize',16),box off,ylim([1 3])
%%
scatterParams.axLim = 'auto'; 
scatterParams.axLen = 6;
scatterParams.plotInExistingFig = 0;
scatterParams.showFanoAll = 0;
scatterParams.mSize = 10;
plotScatter(Result, -100,scatterParams);
text(2.5, 7, '100 ms before target', 'hori', 'center');
plotScatter(Result, 5,scatterParams);
text(2.5, 7, '0 ms before target', 'hori', 'center');
plotScatter(Result, 300, scatterParams);
text(2.5, 7, '100 ms after target', 'hori', 'center');
plotScatter(Result, 600, scatterParams);
text(2.5, 7, '300 ms after target', 'hori', 'center');


%% Plot as mean spike count to variance
figure,
for n = 1:40
    scatter(Result.scatterData(25+n).mn,Result.scatterData(25+n).var,'k','filled'),hold on    
end
for n = 1:40
    scatter(MissResult.scatterData(25+n).mn,MissResult.scatterData(25+n).var,'MarkerEdgeColor','none',...
              'MarkerFaceColor',[.5 .5 .5],...
              'LineWidth',1.5),hold on   
end
line(0:30,0:30,'LineWidth',2,'Color',[ 0 0 0 ])
xlim([0 5])
ylim([0 5])
set(gca,'TickDir','out'),set(gca,'fontsize',12),box off
xlabel('Mean spike rate')
ylabel('Spike variance')
%% Bin data for analysis
[mntotal,valm] = arrayfun(@(x) discretize(x.mn,8),Result.scatterData,'UniformOutput',false);
[vartotal,valvar] = arrayfun(@(x) discretize(x.var,8),Result.scatterData,'UniformOutput',false);
for n = 1:74
    for nn = 1:8
        temp = Result.scatterData(n).mn(mntotal{n}==nn);
        mnSpk(nn,n) = nanmean(temp(temp>0));
        varSpk(nn,n) = nanmean(Result.scatterData(n).var(vartotal{n}==nn));
    end
end
mnSpk = nanmean(mnSpk(:,28:48),2);
varSpk = nanmean(varSpk(:,28:48),2);
varSpk = inpaint_nans(varSpk);
%% 
figure,plot(mnSpk,varSpk)
%%
ScatterMovie(Result);
%%
pad = [];
for n = 20
    Show_Spikes(M1Data(n).spikes);hold on    
end
%% FUNCTIONS
% FFcalculations
function Result = calcNeuralVariance(M1Data)
% Remove empties
temp = find(arrayfun(@(x) isempty(x.spikes), M1Data)==1);
M1Data(temp) = [];
addpath(genpath('C:\Users\khan332\Documents\GitHub\Variance_toolbox'));
% times = 100:15:1200;  % from 200 ms before target onset until 450 ms after.
% fanoParams.alignTime = 500;    % this time will become zero time
% fanoParams.boxWidth = 200;     % 50 ms sliding window.
times = 200:15:1000;  % from 200 ms before target onset until 450 ms after.
fanoParams.alignTime = 500;    % this time will become zero time
fanoParams.boxWidth = 100;     % 50 ms sliding window.
%Result = VarVsMean(M1Data, times, fanoParams);
Result = MeanFano(M1Data, times, fanoParams);
plotFanoParams.plotRawF = 1;
plotFano(Result,plotFanoParams);
end

function FakResult = calcFakNeuralVariace(M1Data)
times = 1000:15:2500;  % from 200 ms before target onset until 450 ms after.
fanoParams.alignTime = 1500;    % this time will become zero time
fanoParams.boxWidth = 100;     % 50 ms sliding window.
FakResult = VarVsMean(Fakerize(M1Data,'poisson'), times, fanoParams);  % takes a while
%FakResult = MeanFano(Fakerize(PMDdata2,'gamma'), times, fanoParams);
plotFanoParams.plotRawF = 1;
plotFano(FakResult, plotFanoParams);
end


function Stats = calcStats(Results)
Var = arrayfun(@(x) horzcat(x.var), Results.scatterData, 'UniformOutput', false);
Var = horzcat(Var{:});
MM = arrayfun(@(x) horzcat(x.mn), Results.scatterData, 'UniformOutput', false);
MM = horzcat(MM{:});
Stats.FF = Var./MM;
Stats.FFdrop = nanmean(Stats.FF(:,1:20),2)-nanmean(Stats.FF(:,21:end),2);
Stats.FFstim = nanmean(Stats.FF(:,21:end),2);
end





