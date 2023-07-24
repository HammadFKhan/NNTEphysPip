function Spikes = spikeRateAnalysis(Spikes)
% Analysis for spike rates across depth which includes basic curation for
% nonresponsive units (ie. units that are essentialy return a 0 for mean
% calculation.

% Input: 
% Spikes - should have a sortedSpike rate function as a function depth
% Output:
% Spikes -  statistics for spike rate as a function of depth

if ~isfield(Spikes,'Depth')
    error('Depth structure does not exist in spikes.')
end
mSpikeRateDepth(:,1) = Spikes.Depth.depth;
mSpikeRateDepth(:,2) = mean(Spikes.Depth.sortedSpikeRate,2);
l23idx = mSpikeRateDepth(:,1)<=300;
l5idx = mSpikeRateDepth(:,1)>300;
L23SpikeRate = mSpikeRateDepth(l23idx,2);
L5SpikeRate = mSpikeRateDepth(l5idx,2);
