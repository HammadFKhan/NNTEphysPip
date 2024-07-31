function Spikes = GPSFAnalysis(Spikes,lfp)
%% Spike field analysis function for lever task
% Here we call the generalized phase spike field function to begin analysis
% across task conditions and cue times of lever task. Along with some basic
% analysis and quantifications for plotting later down the road

% Assert the structures of lfp and spikes are the same. Sometimes the
% windowing parameters may be different or we just loaded the wrong files
% in. 
assert(size(Spikes.PSTH.hit.spks{1},2)==size(squeeze(lfp.hitxgp{1}),2)); %checks length
assert(size(Spikes.PSTH.hit.spks{1},1)==length(lfp.hitxgp));
assert(size(Spikes.PSTH.miss.spks{1},1)==length(lfp.missxgp));
assert(size(Spikes.PSTH.MIHit.spks{1},1)==length(lfp.MIhitxgp));
assert(size(Spikes.PSTH.MIFA.spks{1},1)==length(lfp.MIFAxgp));


% Grab the localized spike electrode so we can match to GP electrode
spkChan = arrayfun(@(x) vertcat(x.channelDepth),Spikes.Clusters)';
precueWin = 1200:1500;
postcueWin = 1501:1800;
% We do pre- post- cue reponse across task function. For hits/miss it'll be
% based on cue response at time 1500. For FA it'll be based on time of
% movement. Which we will also do for motion initiated hits to serve as a
% control for the jittered PSTH
preCuespk = cellfun(@(x) x(:,precueWin), Spikes.PSTH.hit.spks, 'UniformOutput', false);
preCuexgp = cellfun(@(x) x(:,:,precueWin), lfp.hitxgp, 'UniformOutput', false);
postCuespk = cellfun(@(x) x(:,postcueWin), Spikes.PSTH.hit.spks, 'UniformOutput', false);
postCuexgp = cellfun(@(x) x(:,:,postcueWin), lfp.hitxgp, 'UniformOutput', false);


[Spikes.SpikeField.precuehit] = GPSpikeField(preCuespk,preCuexgp,spkChan);
[Spikes.SpikeField.postcuehit] = GPSpikeField(postCuespk,postCuexgp,spkChan);

preCuespk = cellfun(@(x) x(:,precueWin), Spikes.PSTH.miss.spks, 'UniformOutput', false);
preCuexgp = cellfun(@(x) x(:,:,precueWin), lfp.missxgp, 'UniformOutput', false);
postCuespk = cellfun(@(x) x(:,postcueWin), Spikes.PSTH.miss.spks, 'UniformOutput', false);
postCuexgp = cellfun(@(x) x(:,:,postcueWin),lfp.missxgp, 'UniformOutput', false);

[Spikes.SpikeField.precuemiss] = GPSpikeField(preCuespk,preCuexgp,spkChan);
[Spikes.SpikeField.postcuemiss] = GPSpikeField(postCuespk,postCuexgp,spkChan);

preCuespk = cellfun(@(x) x(:,precueWin), Spikes.PSTH.MIFA.spks, 'UniformOutput', false);
preCuexgp = cellfun(@(x) x(:,:,precueWin), lfp.MIFAxgp, 'UniformOutput', false);
postCuespk = cellfun(@(x) x(:,postcueWin), Spikes.PSTH.MIFA.spks, 'UniformOutput', false);
postCuexgp = cellfun(@(x) x(:,:,postcueWin), lfp.MIFAxgp, 'UniformOutput', false);

[Spikes.SpikeField.precueMIFA] = GPSpikeField(preCuespk,preCuexgp,spkChan);
[Spikes.SpikeField.postcueMIFA] = GPSpikeField(postCuespk,postCuexgp,spkChan);

preCuespk = cellfun(@(x) x(:,precueWin), Spikes.PSTH.MIHit.spks, 'UniformOutput', false);
preCuexgp = cellfun(@(x) x(:,:,precueWin), lfp.MIhitxgp, 'UniformOutput', false);
postCuespk = cellfun(@(x) x(:,postcueWin), Spikes.PSTH.MIHit.spks, 'UniformOutput', false);
postCuexgp = cellfun(@(x) x(:,:,postcueWin), lfp.MIhitxgp, 'UniformOutput', false);

[Spikes.SpikeField.precueMIHit] = GPSpikeField(preCuespk,preCuexgp,spkChan);
[Spikes.SpikeField.postcueMIHit] = GPSpikeField(postCuespk,postCuexgp,spkChan);
end
