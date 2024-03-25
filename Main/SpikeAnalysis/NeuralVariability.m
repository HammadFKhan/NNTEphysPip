%% Prep data for VarMean analysis
fpath = 'F:\LeverTask\Ephys\Analysis\spksPooled';
file = dir(fullfile(fpath,'*.mat'));

%%
count = 1;
for fileNum = 1:length(file)
    load(fullfile(file(fileNum).folder,file(fileNum).name))
    for n = 1:length(Spikes.PSTH.hit.spks)
        PMDdata2(count).spikes = logical(Spikes.PSTH.hit.spks{n});
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
temp = find(arrayfun(@(x) isempty(x.spikes), PMDdata2)==1);
PMDdata2(temp) = [];

%%
addpath(genpath('C:\Users\khan332\Documents\GitHub\Variance_toolbox'));
times = 100:25:1200;  % from 200 ms before target onset until 450 ms after.
fanoParams.alignTime = 500;    % this time will become zero time
fanoParams.boxWidth = 80;     % 50 ms sliding window.
%Result = VarVsMean(PMDdata2, times, fanoParams);
Result = MeanFano(PMDdata2, times, fanoParams);
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
plotScatter(Result, 0,scatterParams);
text(2.5, 7, '0 ms before target', 'hori', 'center');
plotScatter(Result, 300, scatterParams);
text(2.5, 7, '100 ms after target', 'hori', 'center');
plotScatter(Result, 600, scatterParams);
text(2.5, 7, '300 ms after target', 'hori', 'center');
%% Check Fakerized data
Result = VarVsMean(Fakerize(PMDdata2,'gamma'), times, fanoParams);  % takes a while
Result = MeanFano(Fakerize(PMDdata2,'gamma'), times, fanoParams);

plotFanoParams.plotRawF = 0;
plotFano(Result, plotFanoParams);
%%
ScatterMovie(Result);
%%
pad = [];
for n = 20
    Show_Spikes(PMDdata2(n).spikes);hold on    
end
