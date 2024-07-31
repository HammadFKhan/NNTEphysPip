function [LaminarData] = GPSpikeField(spikes,lfp,spkChan)
% Now generate structure to get laminar phase
% Code to calculate the laminar profile of linear electrode data. The code
% takes as input:
% spikes: the PSTH times of a single spike across trials. Here we collapse
% the single spike into a Nx1 array to calculate phase
% lfp: is the generalized lfp of complex values in which we calculate phase
% from
% output LaminarData contains the spike-phase distributions, circular
% resultants, and mean phase angles of the spike-LFP relationship on each
% channel.
% Example data
% EXAMPLE INPUT: GPSpikeField(Spikes.PSTH.hit.spks,LFP.probe1.hitxgp,spkChan)
%  10/21/2022 - ZD
%  02/23/2024 - HK
%%% Prep for generalized spike-field coherence
% Here we can use the Spike PSTH and Spike good spike to generate the
% neccessary spike format. 
% This function is called when we run the GPSFAnalysis.m function.

plotOn = 0;
% Double check spk lengths for analysis
assert(length(spkChan)==length(spikes))



NChan = min(size(squeeze(lfp{1})));
NSpikes = length(spkChan);
% Sort spikes from superficial to deep
[~,laminarSpike] = sort(spkChan);
% Initalize Vars
phaseSpikes = cell(NChan,NSpikes);

% Filter and calculate LFP phase
% Chunk data because GP doesn't like big numbers

fprintf('Calculating Spike-field coherence...')
for i = 1:NChan
    phaseLFP = cellfun(@(x) x(i,:), lfp, 'UniformOutput', false);
    phaseLFP = angle(horzcat(phaseLFP{:}));
    for j = 1:NSpikes %For each spike channel
        spikeTimes = reshape(spikes{j},1,[]); % Here we grab the spikes from superficial to deep
        if ~isempty(spikeTimes) %No spikes? No phase.
            theseSpikes = find(spikeTimes==1);
            %Get the spike phases
            phaseSpikes{i,j} = [phaseSpikes{i,j} phaseLFP(1,theseSpikes)];
        else
            phaseSpikes{i,j} = NaN;
        end
    end
end
fprintf('done\n')

fprintf('Calculating preferred spike phase...')
spikePhase = zeros(NChan,NSpikes);
prefPhase = zeros(NChan,NSpikes);
for ii = 1:NChan
    for jj = 1:NSpikes
        spikePhase(ii,jj) = circ_r(phaseSpikes{ii,jj}'); %Calculate circular resultant
        prefPhase(ii,jj) = circ_mean(phaseSpikes{ii,jj}'); %Calculate mean phase angle
    end
end
fprintf('done\n')

fprintf('Calculating lfp-lfp phase...')
phaseCorr = zeros(NChan,NChan);
phaseLFP = cellfun(@squeeze,lfp,'UniformOutput',false);
phaseLFP = angle(horzcat(phaseLFP{:}));
for i = 1:NChan
    for j = 1:NChan
        phaseCorr(i,j) = circ_corrcc(phaseLFP(i,:),phaseLFP(j,:));
    end
end
fprintf('done\n')

%Collect Output
LaminarData.SPI = spikePhase;
LaminarData.Angle = prefPhase;
LaminarData.Dist = phaseSpikes;
LaminarData.PhaseCorr = phaseCorr;
LaminarData.SpikeElectrode = spkChan;

if plotOn
figure
imagesc(spikePhase)
colormap hot
c = colorbar;
c.Label.String = 'SPI';
ylabel('LFP Electrode Channel')
xlabel('Spiking #')
set(gca,'fontsize',14,'linewidth',1.5)

figure
imagesc(prefPhase)
map = colorcet( 'C2' );
map = circshift(map,1);
colormap(map)
c = colorbar;
c.Label.String = 'Best Phase (rad)';
ylabel('LFP Electrode Channel')
xlabel('Spike #')
set(gca,'fontsize',14,'linewidth',1.5)

figure
imagesc(phaseCorr)
colormap hot
c = colorbar;
c.Label.String = 'Circ Corr (r)';
ylabel('Channel')
xlabel('Channel')
set(gca,'fontsize',14,'linewidth',1.5)
end