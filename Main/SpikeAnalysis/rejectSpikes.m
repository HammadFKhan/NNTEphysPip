function Spikes = rejectSpikes(Spikes,fractionTrials,cutoffFR,parameters)
% Function to reject all the spikes that do not satisfy the following
% criterion 
% 1. Spike fires atleast in 75% (fractionTrials) of the HIT trials (Should I include Hits and FAs?)
% 2. Either pre-stimulus or post-stimulus firing rate is more than 1Hz (cutoffFR)
% 3. Example input: Spikes = rejectSpikes(Spikes,0.75,1,parameters)

Spikes.rawClusters = Spikes.Clusters;
Spikes.rawPSTH = Spikes.PSTH;

%% Rejecting spikes that do not fire for fractionTrials % of HIT trials (Should I include Hits and FAs?)
Spikes.badSpikes = [];

for i=1:size(Spikes.PSTH.hit.spks,2) % interating over number of spikes
    nSpikesPerTrial = sum(Spikes.PSTH.hit.spks{1,i},2);
    spikePerTrial = find(nSpikesPerTrial);
    if numel(spikePerTrial)<fractionTrials*size(Spikes.rawPSTH.hit.spks{1,1},1)
        Spikes.badSpikes = [Spikes.badSpikes i];
    end
end

disp(['Rejected ' num2str(numel(Spikes.badSpikes)) ' Spike Clusters based on fraction criterion']);
Spikes.Clusters(Spikes.badSpikes) = [];
Spikes.PSTH.hit.spks(Spikes.badSpikes) = [];
Spikes.PSTH.hit.spkRates(Spikes.badSpikes,:) = [];
Spikes.PSTH.miss.spks(Spikes.badSpikes) = [];
Spikes.PSTH.miss.spkRates(Spikes.badSpikes,:) = [];
Spikes.PSTH.MIHit.spks(Spikes.badSpikes) = [];
Spikes.PSTH.MIHit.spkRates(Spikes.badSpikes,:) = [];
Spikes.PSTH.MIFA.spks(Spikes.badSpikes) = [];
Spikes.PSTH.MIFA.spkRates(Spikes.badSpikes,:) = [];

%% Rejecting spikes that do not have FR > cutoffFR either pre or post stimulus 
% (JUST FOR HIT TRIALS - Logic it will have both movement and cue responsive neurons)
Spikes.badSpikes2 = [];

for i=1:size(Spikes.PSTH.hit.spks,2) % interating over number of spikes
    preFR = (sum(Spikes.PSTH.hit.spks{1,i}(:,1:parameters.windowBeforeCue*parameters.Fs),"all")/size(Spikes.rawPSTH.hit.spks{1,1},1))/parameters.windowBeforeCue;
    postFR = (sum(Spikes.PSTH.hit.spks{1,i}(:,parameters.windowBeforeCue*parameters.Fs+1:(parameters.windowBeforeCue+parameters.windowAfterCue)*parameters.Fs+1),"all")/size(Spikes.rawPSTH.hit.spks{1,1},1))/parameters.windowAfterCue;
    if (preFR < cutoffFR || postFR < cutoffFR)
        Spikes.badSpikes2 = [Spikes.badSpikes2 i];
    end
end

disp(['Rejected ' num2str(numel(Spikes.badSpikes2)) ' Spikes Clusters based on FR critereon']);
Spikes.Clusters(Spikes.badSpikes2) = [];
Spikes.PSTH.hit.spks(Spikes.badSpikes2) = [];
Spikes.PSTH.hit.spkRates(Spikes.badSpikes2,:) = [];
Spikes.PSTH.miss.spks(Spikes.badSpikes2) = [];
Spikes.PSTH.miss.spkRates(Spikes.badSpikes2,:) = [];
Spikes.PSTH.MIHit.spks(Spikes.badSpikes2) = [];
Spikes.PSTH.MIHit.spkRates(Spikes.badSpikes2,:) = [];
Spikes.PSTH.MIFA.spks(Spikes.badSpikes2) = [];
Spikes.PSTH.MIFA.spkRates(Spikes.badSpikes2,:) = [];

%% Sorting good spikes according to depth
[~,sortspikeDepth] = sort(cell2mat(arrayfun(@(s) s.channelDepth, Spikes.Clusters,'UniformOutput',false))');
Spikes.Clusters = Spikes.Clusters(sortspikeDepth);
Spikes.nSpikes = size(Spikes.Clusters,2);

disp([num2str(Spikes.nSpikes) ' number of good spike clusters detected']);