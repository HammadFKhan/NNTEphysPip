function prepWrap(Spikes,ds_filename)
% Prep data for wrapping structure to load into python. 
SqSpikes = zeros(size(Spikes.PSTH.hit.spks{1},1),size(Spikes.PSTH.hit.spks{1},2),length(Spikes.PSTH.hit.spks));
for n = 1:size(SqSpikes,3)
    SqSpikes(:,:,n) = Spikes.PSTH.hit.spks{n};
end
pullIndex = Spikes.PSTH.hit.pl;
pull1 = pullIndex(:,1);
pull2 = pullIndex(:,2);
pull3 = pullIndex(:,3);
tmin = 0;
tmax = size(SqSpikes,2)-1;

% Trial IDs, spiketimes, and neuron_ids are a flatten matrix of the 3d
% array. Trial IDs and spike times should have zero indexing
trial_ids = [];
spiketimes = [];
neuron_ids = [];
for ntrials = 1:size(SqSpikes,1)
    spkTemp = squeeze(SqSpikes(ntrials,:,:))';
    [r,c,v] = find(spkTemp==1);
    trial_ids = [trial_ids;(ntrials-1)*ones(size(r))];
    spiketimes = [spiketimes;c];
    neuron_ids = [neuron_ids;r];
end
neuron_ids = neuron_ids-1; % zero indexing
spiketimes = spiketimes-1; % zero indexing
[fpath,name,exts] = fileparts(ds_filename);
sessionName = [fpath,'\','spikes_to_Warp.mat'];
save(sessionName,"tmin","tmax","pull1","pull2","pull3","trial_ids","spiketimes","neuron_ids","fpath");
disp('Data Saved for warping')