%%
clear
%fpath = 'F:\LeverTask\Ephys\Analysis\M2Spikes';
%fpath = 'F:\LeverTask\Ephys\Analysis\spksPooledwFA'
fpath = 'D:\M1_GSP';
%fpath = 'D:\M2SpikeData';
file = dir(fullfile(fpath,'*.mat'));

%%%
M2Datahit = struct();
M2Datamiss = struct();
M2DataMIFA = struct();
count1 = 1;count2 = 1;count3 = 1;

totalSpikes.hit = [];
totalSpikes.miss = [];
totalSpikes.FA = [];
precisionHitSpk = nan(1,length(file));
precisionMissSpk = nan(1,length(file));
precisionFASpk = nan(1,length(file));
SqEntropyMice = struct();
for fileNum = 1:length(file)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(file(fileNum).folder,file(fileNum).name))
    totalSpikes.hit = [totalSpikes.hit;Spikes.PSTH.hit.normSpk];
    totalSpikes.miss = [totalSpikes.miss;Spikes.PSTH.miss.normSpk];
    totalSpikes.FA = [totalSpikes.FA;Spikes.PSTH.MIFA.normSpk];
    
    dat = Spikes.PSTH.hit.normSpk;
    precisionHitSpk(fileNum) = nanmean(nanmean(calcSpkPrecision(dat,0),2));
    dat = Spikes.PSTH.miss.normSpk;
    precisionMissSpk(fileNum) = nanmean(nanmean(calcSpkPrecision(dat,0),2));
    dat = Spikes.PSTH.MIFA.normSpk;
    precisionFASpk(fileNum) = nanmean(nanmean(calcSpkPrecision(dat,0),2));
    %SqEntropyMice(fileNum).SqEntropy = getSqEntropy(Spikes);
end
%%
normSpikeRate = totalSpikes.hit;
f = figure,subplot(131)
plotSpkSeq(normSpikeRate)
%title('Hit')
colormap(flip(gray))

set(gca,'fontsize',16)
caxis([0.0 2])
xlim([-500 1500])
ylim([0 525])
normSpikeRate = totalSpikes.miss;
subplot(132),plotSpkSeq(normSpikeRate)
%title('Miss')
colormap(flip(gray))
xlim([-500 1500])
ylim([0 525])
set(gca,'fontsize',16)
caxis([0.0 2])
normSpikeRate = totalSpikes.FA;
subplot(133),plotSpkSeq(normSpikeRate)
%title('FA')
colormap(flip(gray))
set(gca,'fontsize',16)
caxis([0.0 2])
xlim([-500 1500])
ylim([0 525])
f.Position = [681 559 560 250];
%% Calculate precision of sequences
% Taken from: Changes in the neural control of a complex motor sequence during learning
% Bence P. Ölveczky,Timothy M. Otchy,Jesse H. Goldberg,Dmitriy Aronov, and Michale S. Fee
% https://doi.org/10.1152/jn.00018.2011
% 
% The precision of the song-aligned spike trains was measured using average
% pairwise correlation across all pairs of spike train for a given
% condition. Spike trains were converted into instantaneous firing rates
% R(t) as follows:
% 
% 𝑅⁡(𝑡)= 1 𝑡𝑖+1−𝑡𝑖
%  
% ;	for	𝑡𝑖<𝑡≤𝑡𝑖+1
%  ,
% where ti is the ith spike. These instantaneous firing rates were then
% convolved with a 8-ms Gaussian function (Leonardo and Fee 2005), yielding
% a smoothed firing rate function r(t).
% The correlation coefficient (CC) was then calculated between these firing
% rate functions for all pairs of spike trains as follows: CC= 1 𝑁pairs
%  
% ⁢ 𝑁 ∑ 𝑖 𝑁 ∑ 𝑗>𝑖 CC𝑖⁢𝑗, CC𝑖⁢𝑗= ⟨ ̂ 𝑟 𝑖⁡(𝑡)⋅ ̂ 𝑟 𝑗⁡(𝑡)⟩𝑡 √⟨
% ̂ 𝑟 𝑖⁡(𝑡)2⟩𝑡⁢⟨ ̂ 𝑟 𝑗⁡(𝑡)2⟩𝑡
%  
% , where r̂(t) is the mean-subtracted smoothed firing rate function.

%%% CALCULATE ACROSS SESSIONS
dat = [precisionHitSpk',precisionMissSpk',precisionFASpk'];
figure,customBarplot(dat,'Scatter','off'),hold on

for n = 1:size(dat,1)
    scatter(1:3,[dat(n,1),dat(n,2),dat(n,3)],'k','filled')
    line([1 2 3],[dat(n,1),dat(n,2),dat(n,3)])
end

ylabel('Spike Precision')
box off, set(gca,'tickdir','out','fontsize',14)
[~,~,stats] = anova1([precisionHitSpk',precisionMissSpk',precisionFASpk'])
c = multcompare(stats)



%%
%%% FUNCTION CALL

function plotSpkSeq(normSpikeRate)
idx = zeros(size(normSpikeRate,1),1);
for n = 1:length(idx)
    [~,idx(n)] = max(normSpikeRate(n,:));
end
[~,idxc] = sort(idx);
%%%
path = idx(idxc);
imagesc(-1.5*1000:1.5*1000,...
    1:size(normSpikeRate,1),normSpikeRate(idxc,:)),hold on
plot((path)-1.5*1000,1:size(normSpikeRate,1),'r','LineWidth',1)
end

function precisionSpk = calcSpkPrecision(dat,plotOn)
precisionSpk = abs(corrcoef(dat'));
if plotOn
    figure,imagesc(precisionSpk),caxis([0 1]),axis square, colorbar
    precisionSpk(precisionSpk==1) = NaN;
    figure,histogram(precisionSpk,0:0.01:1,'Normalization','Probability','EdgeColor','none')
    box off, set(gca,'tickdir','out','fontsize',14),axis square
    xlabel('Spike Precision'),ylabel('Probability')
end

end