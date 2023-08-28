%% test custom automerge
% test scipt for automergers, checks if spike can be merged based on
% similarity index. Can also further refine response by parsing if unit has
% been labeled as good isolated unit or MUA by kilosort
simIndx = readNPY(fullfile(path, 'similar_templates.npy'));
%%
% Remove 0 indexing for brevity
SpikeClusterssample = SpikeClusters+1;
% Go across dimensions (technically half efficient since we'll override
% previous if sim is same)
for n = 1:length(simIndx)
    mergeCan = (find(simIndx(n,:)>0.8));
    for nn = 1:length(mergeCan) % Parse through indexs that meet criteria
        SpikeClusterssample(SpikeClusterssample==mergeCan(nn)) = n; 
    end
end
% rewrite as new spikeClusters
writeNPY(uint32(SpikeClusterssample-1), fullfile(path, 'spike_clusters.npy')); % -1 for zero indexing
