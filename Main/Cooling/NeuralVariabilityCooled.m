%% Prep data for VarMean analysis
%fpath = 'F:\LeverTask\Ephys\Analysis\spksPooled';
%fpath = 'F:\LeverTask\Ephys\Analysis\spksPooledwFA';
%fpath = 'F:\LeverTask\Ephys\Analysis\M2Spikes';
fpath = 'F:\LeverTask\Ephys\Cooling\Spikes';
file = dir(fullfile(fpath,'*.mat'));

%%%
M1Datahit = struct();
M1DatahitCooled = struct();
M1Datamiss = struct();
M1DatamissCooled = struct();
M1DataMIFA = struct();
M1DataMIFACooled = struct();
count1 = 1;count2 = 1;count3 = 1;


for fileNum = 1:length(file)
    load(fullfile(file(fileNum).folder,file(fileNum).name))
    
    for n = 1:IntanBehaviour.nCueHit
        IntanBehaviour.hitTemp(n,1) = IntanBehaviour.temperature(IntanBehaviour.cueHitTrace(n).LFPIndex(1));
    end
    IntanBehaviour.hitTemp = IntanBehaviour.hitTemp-IntanBehaviour.temperature(100);
    for n = 1:IntanBehaviour.nCueMiss
        IntanBehaviour.missTemp(n,1) = IntanBehaviour.temperature(IntanBehaviour.cueMissTrace(n).LFPIndex(1));
    end
    IntanBehaviour.missTemp = IntanBehaviour.missTemp-IntanBehaviour.temperature(100);
    for n = 1:length(IntanBehaviour.missTrace)
        IntanBehaviour.FATemp(n,1) = IntanBehaviour.temperature(IntanBehaviour.missTrace(n).LFPIndex(1));
    end
    IntanBehaviour.FATemp = IntanBehaviour.FATemp-IntanBehaviour.temperature(100);
    
    coolHitTrials = IntanBehaviour.hitTemp<-10;
    coolMissTrials = IntanBehaviour.missTemp<-10;
    coolFATrials = IntanBehaviour.FATemp<-10;
    %%% Segment as a function of temp. Here we consider all trial
    %%% conditions were probe temperature was -10 C from baseline.
    for n = 1:length(Spikes.PSTH.hit.spks)
        M1Datahit(count1).spikes = logical(Spikes.PSTH.hit.spks{n}(~coolHitTrials,:));
        M1DatahitCooled(count1).spikes = logical(Spikes.PSTH.hit.spks{n}(coolHitTrials,:));
        count1 = count1+1;
    end
    for n = 1:length(Spikes.PSTH.miss.spks)
        M1Datamiss(count2).spikes = logical(Spikes.PSTH.miss.spks{n}(~coolMissTrials,:));
        M1DatamissCooled(count2).spikes = logical(Spikes.PSTH.miss.spks{n}(coolMissTrials,:));
        count2 = count2+1;
    end
    for n = 1:length(Spikes.PSTH.MIFA.spks)
        M1DataMIFA(count3).spikes = logical(Spikes.PSTH.MIFA.spks{n}(~coolFATrials,:));
        M1DataMIFACooled(count3).spikes = logical(Spikes.PSTH.MIFA.spks{n}(coolFATrials,:));
        count3 = count3+1;
    end
end
%% Calculate Neural variability via mean-matched FF
% See: Churchland MM et al. (2010)  Stimulus onset quenches neural variability: a widespread cortical phenomenon. Nat. Neurosci.

% DEFINE Data type here
% Call FF function and calculate
M1Data = M1Datahit;
Results.Baseline.Hit = calcNeuralVariance(M1Data);
%%
M1Data = M1Datamiss;
Results.Baseline.miss = calcNeuralVariance(M1Data);
M1Data = M1DataMIFA;
Results.Baseline.FA = calcNeuralVariance(M1Data);
%%
% Call FF function and calculate but for cooled data
M1Data = M1DatahitCooled;
Results.Cooled.Hit = calcNeuralVariance(M1Data);
%%
M1Data = M1DatamissCooled;
Results.Cooled.miss = calcNeuralVariance(M1Data);
M1Data = M1DataMIFACooled;
Results.Cooled.FA = calcNeuralVariance(M1Data);
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

M1Data = M1Datahit;
Results.Fak.hit = calcFakNeuralVariace(M1Data);
M1Data = M1Datamiss;
Results.Fak.miss = calcFakNeuralVariace(M1Data);
M1Data = M1DataMIFA;
Results.Fak.FA = calcFakNeuralVariace(M1Data);
%% Calc FF stats and make structure for further analysis
Stats.Baseline.Hit = calcStats(Results.Baseline.Hit);
Stats.Baseline.Miss = calcStats(Results.Baseline.miss);
Stats.Baseline.FA = calcStats(Results.Baseline.FA);

Stats.Cooled.Hit = calcStats(Results.Cooled.Hit);
Stats.Cooled.Miss = calcStats(Results.Cooled.miss);
Stats.Cooled.FA = calcStats(Results.Cooled.FA);

Stats.Fak.Hit = calcStats(Results.Fak.hit);
Stats.Fak.Miss = calcStats(Results.Fak.miss);
Stats.Fak.FA = calcStats(Results.Fak.FA);
%% Plot it out
figure(1)
clf
dat = Stats.Baseline;
Stats.Baseline = plotFFStats(dat);
figure(2)
clf
dat = Stats.Cooled;
Stats.Cooled = plotFFStats(dat);
figure(3)
clf
dat = Stats.Fak;
Stats.Fak = plotFFStats(dat);
%% Stats
[~,~,Stats.Baseline.statsFFdrop] = anova1(Stats.Baseline.FFdrop);
Stats.Baseline.resultsFFdrop = multcompare(Stats.Baseline.statsFFdrop);
[~,~,Stats.Baseline.statsFFstim] = anova1(Stats.Baseline.FFstim);
Stats.Baseline.resultsFFstim = multcompare(Stats.Baseline.statsFFstim);

[~,~,Stats.Cooled.statsFFdrop] = anova1(Stats.Cooled.FFdrop);
Stats.Cooled.resultsFFdrop = multcompare(Stats.Cooled.statsFFdrop);
[~,~,Stats.Cooled.statsFFstim] = anova1(Stats.Cooled.FFstim);
Stats.Cooled.resultsFFstim = multcompare(Stats.Cooled.statsFFstim);

[~,~,Stats.Fak.statsFFdrop] = anova1(Stats.Fak.FFdrop);
Stats.Fak.resultsFFdrop = multcompare(Stats.Fak.statsFFdrop);
[~,~,Stats.Fak.statsFFstim] = anova1(Stats.Fak.FFstim);
Stats.Fak.resultsFFstim = multcompare(Stats.Fak.statsFFstim);
%%
[~,~,stats] = anova1(dat);
multcompare(stats);
%% Save Data
sessionName = [fpath,'/','Stats.mat'];
% save(sessionName,"IntanBehaviour","fpath","parameters","-v7.3");
save(sessionName,"Stats","Results","fpath","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
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
% fanoParams.boxWidth = 100;     % 50 ms sliding window.
times = 1000:15:2500;  % from 200 ms before target onset until 450 ms after.
fanoParams.alignTime = 1500;    % this time will become zero time
fanoParams.boxWidth = 100;     % 50 ms sliding window.
%Result = VarVsMean(M1Data, times, fanoParams);
Result = MeanFano(M1Data, times, fanoParams);
plotFanoParams.plotRawF = 1;
plotFano(Result,plotFanoParams);
end

function FakResult = calcFakNeuralVariace(M1Data)
% times = 100:15:1200;  % from 200 ms before target onset until 450 ms after.
% fanoParams.alignTime = 500;    % this time will become zero time
% fanoParams.boxWidth = 200;     % 50 ms sliding window.

times = 1000:15:2500;  % from 200 ms before target onset until 450 ms after.
fanoParams.alignTime = 1500;    % this time will become zero time
fanoParams.boxWidth = 100;     % 50 ms sliding window.
%FakResult = VarVsMean(Fakerize(M1Data,'poisson'), times, fanoParams);  % takes a while
FakResult = MeanFano(Fakerize(M1Data,'gamma'), times, fanoParams);
plotFanoParams.plotRawF = 1;
plotFano(FakResult, plotFanoParams);
end


function Stats = calcStats(Results)
Var = arrayfun(@(x) horzcat(x.var), Results.scatterData, 'UniformOutput', false);
Var = horzcat(Var{:});
MM = arrayfun(@(x) horzcat(x.mn), Results.scatterData, 'UniformOutput', false);
MM = horzcat(MM{:});
Stats.FF = Var./MM;
Stats.FFdrop = nanmean(Stats.FF(:,1:34),2)-nanmean(Stats.FF(:,35:70),2);
Stats.FFstim = nanmean(Stats.FF(:,35:70),2);
end

function output = plotFFStats(dat)

colors = [0 0.4470 0.7410;0.75 0.75 0.75;190/255 30/255 45/255];

dat1 = dat.Hit.FFdrop;
dat2 = dat.Miss.FFdrop;
dat3 = dat.FA.FFdrop;

temp = nan(max([length(dat1) length(dat2) length(dat3)]),3);
temp(1:length(dat1),1) = dat1;
temp(1:length(dat2),2) = dat2;
temp(1:length(dat3),3) = dat3;
subplot(121),violinplot(temp,[],'ShowData',true,'ShowWhiskers',false,'ShowBox',false,'MarkerSize',5,'ViolinColor',colors);
set(gca,'tickdir','out','fontsize',16),box off,ylim([-1.5 1.5]),ylabel('post-cue fanofactor')
output.FFdrop = temp;
dat4 = dat.Hit.FFstim;
dat5 = dat.Miss.FFstim;
dat6 = dat.FA.FFstim;

temp = nan(max([length(dat4) length(dat5) length(dat6)]),3);
temp(1:length(dat4),1) = dat4;
temp(1:length(dat5),2) = dat5;
temp(1:length(dat6),3) = dat6;
subplot(122),violinplot(temp,[],'ShowData',true,'ShowWhiskers',false,'ShowBox',false,'MarkerSize',5,'ViolinColor',colors);
set(gca,'tickdir','out','fontsize',16),box off,ylim([1 3]),ylabel('cue fanofactor')
output.FFstim = temp;
end