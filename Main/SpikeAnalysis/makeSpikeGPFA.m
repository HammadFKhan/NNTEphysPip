function Spikes = makeSpikeGPFA(Spikes, rebalanceFlag)
% formats kilsort spike structure for GPFA analysis each spike is binned
% across 1ms for each column. A new structure row is made for each trial
% The spikes generated from leverPSTH.m is used in this procedure
if ~exist('rebalanceFlag','var')
    rebalanceFlag = 0;
end
if rebalanceFlag
    disp('Rebalanced PSTH flag has been set!')
    Spikes = rebalanceSpikes(Spikes);
    Spikes.unbalancedPSTH = Spikes.PSTH;
    Spikes.PSTH = Spikes.BalancedPSTH;
end
if ~isfield(Spikes,'PSTH'), error('Run leverPSTH.m first!');end

if isfield(Spikes.PSTH,'hit')
    output = organizeSpikes(Spikes.PSTH.hit.spks);
else
    output = [];
end
Spikes.GPFA.hit.dat = output;
if isfield(Spikes.PSTH,'miss')
    output = organizeSpikes(Spikes.PSTH.miss.spks);
else
    output = [];
end
Spikes.GPFA.miss.dat = output;

if isfield(Spikes.PSTH,'MIHit')
    output = organizeSpikes(Spikes.PSTH.MIHit.spks);
    Spikes.GPFA.MIHit.dat = output;
else
    output = [];
    Spikes.GPFA.MIHit.dat = output;
end
if isfield(Spikes.PSTH,'MIFA')
    output = organizeSpikes(Spikes.PSTH.MIFA.spks);
    Spikes.GPFA.MIFA.dat = output;
end


if isfield(Spikes.PSTH,'effortperturb')
    output = organizeSpikes(Spikes.PSTH.effortperturb.spks);
    Spikes.GPFA.effortperturb.dat = output;
end


% if isfield(Spikes.BiPOLES,'hit')
%     output = organizeSpikes(Spikes.BiPOLES.hit.spks);
%     Spikes.GPFA.hit.dat = output;
%     warning("Hit spikes overwritten by BiPOLES")
% end

if rebalanceFlag
    Spikes.PSTH = Spikes.unbalancedPSTH;
    Spikes = rmfield(Spikes,'unbalancedPSTH');
end
%% ------------functions---------
%% Outcomes you want to balance
    function Spikes = rebalanceSpikes(Spikes)
        outcomes = {'hit','miss','MIHit','MIFA'};

        % Assume same neuron count and trial ordering across outcome PSTHs
        nNeurons = numel(Spikes.PSTH.hit.spks);

        %%% 1. Find trial counts for each outcome
        nTrialsPerOutcome = zeros(numel(outcomes),1);
        for oi = 1:numel(outcomes)
            thisField = outcomes{oi};
            % one neuron (e.g. 1) is enough to get nTrials, assuming consistent
            tempSpk = Spikes.PSTH.(thisField).spks{1};   % [nTrials x nTime]
            nTrialsPerOutcome(oi) = size(tempSpk,1);
        end

        %%% 2. Choose common balanced trial count
        % simplest: use the minimum across outcomes
        balTrials = min(nTrialsPerOutcome);

        % If you want to cap at a smaller value (e.g. 40):
        % balTrials = min(40, min(nTrialsPerOutcome));

        %%% 3. Draw balanced trial indices for each outcome
        % Store indices for later reuse across neurons
        balIdx = cell(numel(outcomes),1);
        rng(0);  % reproducible

        for oi = 1:numel(outcomes)
            nThis   = nTrialsPerOutcome(oi);
            % random subset of trials for this outcome
            balIdx{oi} = randsample(nThis, balTrials, false);  % [balTrials x 1]
        end

        %%% 4. Build a new balanced PSTH structure (same layout as Spikes.PSTH)
        Spikes.BalancedPSTH = struct();  % new container

        for oi = 1:numel(outcomes)
            thisField = outcomes{oi};
            thisIdx   = balIdx{oi};              % trial indices to keep
            balSpk    = cell(1, nNeurons);

            for neuronId = 1:nNeurons
                tempSpk = Spikes.PSTH.(thisField).spks{neuronId};  % [nTrials x nTime]
                if isempty(tempSpk)
                    balSpk{neuronId} = [];
                else
                    balSpk{neuronId} = tempSpk(thisIdx, :);        % [balTrials x nTime]
                end
            end

            % copy into new balanced structure
            Spikes.BalancedPSTH.(thisField).spks = balSpk;
        end
    end

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