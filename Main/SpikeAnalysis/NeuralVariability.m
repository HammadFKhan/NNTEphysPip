%% Prep data for VarMean analysis
%fpath = 'F:\LeverTask\Ephys\Analysis\spksPooledwFA';
fpath = 'F:\LeverTask\Ephys\Cooling\Spikes';
file = dir(fullfile(fpath,'*.mat'));

%%%
M1Data = struct();
count = 1;
for fileNum = 1:length(file)
    load(fullfile(file(fileNum).folder,file(fileNum).name))
    for n = 1:length(Spikes.PSTH.hit.spks)
        M1Data(count).spikes = logical(Spikes.PSTH.hit.spks{n});
        count = count+1;
    end
end
% for nn = 1:length(Spikes.PSTH.miss.spks)
%     PMDdata2(n+nn).spikes = logical(Spikes.PSTH.miss.spks{n});
% end
% for nnn = 1:length(Spikes.PSTH.MIFA.spks)
%     PMDdata2(n+nn+nnn).spikes = logical(Spikes.PSTH.MIFA.spks{n});
% end
% Remove empties
temp = find(arrayfun(@(x) isempty(x.spikes), M1Data)==1);
M1Data(temp) = [];

%%
addpath(genpath('C:\Users\khan332\Documents\GitHub\Variance_toolbox'));
% times = 100:15:1200;  % from 200 ms before target onset until 450 ms after.
% fanoParams.alignTime = 500;    % this time will become zero time
% fanoParams.boxWidth = 200;     % 50 ms sliding window.
times = 1000:15:2500;  % from 200 ms before target onset until 450 ms after.
fanoParams.alignTime = 1500;    % this time will become zero time
fanoParams.boxWidth = 200;     % 50 ms sliding window.
Result = VarVsMean(M1Data, times, fanoParams);
%Result = MeanFano(M1Data, times, fanoParams);
plotFanoParams.plotRawF = 1;
plotFano(Result,plotFanoParams);
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

%% Check Fakerized data
FakResult = VarVsMean(Fakerize(M1Data,'gamma'), times, fanoParams);  % takes a while
%FakResult = MeanFano(Fakerize(PMDdata2,'gamma'), times, fanoParams);

plotFanoParams.plotRawF = 1;
plotFano(FakResult, plotFanoParams);
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
