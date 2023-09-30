%% combine spike maps
[path] = uigetdir();
fs  = dir(fullfile(path,'*.mat'));
%%
hl23spk = [];ml23spk = [];hl5spk = [];ml5spk = [];
for n = 1:length(fs)
load(fullfile(fs(n).folder,fs(n).name));
[hitl23,missl23,hitl5,missl5] = sortspkLayer(Spikes);
hl23spk = vertcat(hl23spk,hitl23);
ml23spk = vertcat(ml23spk,missl23);
hl5spk = vertcat(hl5spk,hitl5);
ml5spk = vertcat(ml5spk,missl5);
end
%%
l23Spikes.PSTH.hit.spkRates = hl23spk;
l23Spikes.PSTH.miss.spkRates = ml23spk;
l5Spikes.PSTH.hit.spkRates = hl5spk;
l5Spikes.PSTH.miss.spkRates = ml5spk;
%%
[l23Spikes] = sortSpkLever(l23Spikes,IntanBehaviour);
l5Spikes = sortSpkLever(l5Spikes,IntanBehaviour);
%%
temp = totalSpikes.PSTH.hit.normSpk;
figure,imagesc(temp),colormap(flip(gray)),caxis([0 2.56])
%%%
function [hitl23spks,missl23spks,hitl5spks,missl5spks] = sortspkLayer(Spikes)
l23idx = arrayfun(@(x) x.spikeDepth<350,Spikes.Clusters,'UniformOutput',false); % L23
l5idx = arrayfun(@(x) x.spikeDepth>=350,Spikes.Clusters,'UniformOutput',false); % L5
l23idx = find(cell2mat(l23idx)==1);
l5idx = find(cell2mat(l5idx)==1);
hitl23spks = Spikes.PSTH.hit.spkRates(l23idx,:);
hitl5spks = Spikes.PSTH.hit.spkRates(l5idx,:);
missl23spks = Spikes.PSTH.miss.spkRates(l23idx,:);
missl5spks = Spikes.PSTH.miss.spkRates(l5idx,:);
end

function output = findLayer(spkIdx,idx)
output = ismember(spkIdx,idx);
output = spkIdx(output);
end