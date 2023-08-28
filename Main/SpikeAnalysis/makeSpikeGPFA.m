function Spikes = makeSpikeGPFA(Spikes)
% formats kilsort spike structure for GPFA analysis each spike is binned
% across 1ms for each column. A new structure row is made for each trial
% The spikes generated from leverPSTH.m is used in this procedure

if ~isfield(Spikes,'PSTH'), error('Run leverPSTH.m first!');end

if isfield(Spikes.PSTH,'hit')
    output = organizeSpikes(Spikes.PSTH.hit);
end
Spikes.GPFA.hit.dat = output;
if isfield(Spikes.PSTH,'miss')
    output = organizeSpikes(Spikes.PSTH.miss);
end
Spikes.GPFA.miss.dat = output;
%% ------------functions---------
    function output = organizeSpikes(spk)
        output = struct();
        trials = size(spk{1},1);
        for n = 1:trials
            spikePop = cellfun(@(x) x(n,:),spk,'UniformOutput',false);
            output(n).trialId = n;
            output(n).spikes = vertcat(spikePop{:});
        end
    end
end