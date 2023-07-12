function Spikes = spikeRateAnalysis(Spikes)
% Analysis for spike rates across depth which includes basic curation for
% nonresponsive units (ie. units that are essentialy return a 0 for mean
% calculation.

% Input: 
% Spikes - should have a sortedSpike rate function as a function depth
% Output:
% Spikes -  statistics for spike rate as a function of depth

if ~isfield(Spikes,'spikeRate')
    error('Depth structure does not exist in spikes.')
end
L23SR = Spikes.spikeRate.L23SR;
L23SR(L23SR==0) = NaN;
L4SR = Spikes.spikeRate.L4SR;
L4SR(L4SR==0) = NaN;
L5SR = Spikes.spikeRate.L5SR;
L5SR(L5SR==0) = NaN;
L5ABSR = Spikes.spikeRate.L5ABSR;
L5ABSR(L5ABSR==0) = NaN;
Spikes.spikeRate.mL23SR = nanmedian(L23SR,2);
Spikes.spikeRate.mL4SR = nanmedian(L4SR,2);
Spikes.spikeRate.mL5SR = nanmedian(L5SR,2);
Spikes.spikeRate.mL5ABSR = nanmedian(L5ABSR,2);
end