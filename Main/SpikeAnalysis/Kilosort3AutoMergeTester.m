%% test custom automerge
% test scipt for automergers, checks if spike can be merged based on
% similarity index. Can also further refine response by parsing if unit has
% been labeled as good isolated unit or MUA by kilosort
spikeClusters = readNPY(fullfile(path,'spike_clusters.npy'));
Xsim = readNPY(fullfile(path, 'similar_templates.npy'));
%
% Remove 0 indexing for brevity
SpikeClusterssample = spikeClusters+1;
% Go across dimensions (technically half efficient since we'll override
% previous if sim is same)
Xsim = Xsim - diag(diag(Xsim));
for n = 1:length(Xsim)
    mergeCan = (find(Xsim(n,:)>0.85)); %%% 85% similarity
    for nn = 1:length(mergeCan) % Parse through indexs that meet criteria
        if mergeCan(nn)>n
            SpikeClusterssample(SpikeClusterssample==mergeCan(nn)) = n;
        else
            break;
        end
    end
end
% now apply adjusted cluster properties to variables
mergeId = double(unique(SpikeClusterssample)); % identify new ids
% Include merged templates
% temps = temps(ismember(initClust,mergeId),:,:);
% tempsind = tempsind(ismember(initClust,mergeId),:,:);
% simIndx = simIndx(ismember(initClust,mergeId),ismember(initClust,mergeId));
% make new folder and rewrite spike files for phy2
if ~exist([path '/merged'],'dir')
    mkdir([path '/merged']);
end
writeNPY(uint32(SpikeClusterssample-1), fullfile([path '\merged'], 'spike_clusters.npy')); % -1 for zero indexing
writeNPY(uint32(SpikeClusterssample-1), fullfile([path '\merged'], 'spike_templates.npy')); % -1 for zero indexing
writeNPY(single(Xsim), fullfile([path '\merged'], 'similar_templates.npy'));
copyfile([path,'/params.py'],[path,'/merged']);
copyfile([path,'/amplitudes.npy'],[path,'/merged']);
copyfile([path,'/channel_map.npy'],[path,'/merged']);
copyfile([path,'/channel_positions.npy'],[path,'/merged']);
copyfile([path,'/spike_times.npy'],[path,'/merged']);
copyfile([path,'/templates.npy'],[path,'/merged']);
copyfile([path,'/templates_ind.npy'],[path,'/merged']);
copyfile([path,'/whitening_mat.npy'],[path,'/merged']);
copyfile([path,'/whitening_mat_inv.npy'],[path,'/merged']);






