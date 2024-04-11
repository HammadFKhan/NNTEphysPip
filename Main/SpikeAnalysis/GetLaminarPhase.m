function [LaminarData] = GetLaminarPhase(spikeTimes,LFP,chanMap,Fs)
% Code to calculate the laminar profile of linear electrode data. The code
% takes as input:
% spikeTimes = a Mx1 cell array of M channels having N spike times in msec
% in each cell.
% LFP = a MxN array of M channels having N time points of LFP.
% chanMap is the decending order of channel positions from the surface
% If the data is already in it's physical order, just input an array
% from 1:N Channels
% Fs = the sampling rate per second (I use LFP @ 1000Hz)
% filtBand = the lower and upper bands for the bandpass filter (the actual
% filter settings don't seem to matter too much, but lower and wider seem better)
% GP = Use generalized phase? 1 = yes, 0 = no, use Hilbert (not a big
% difference unless broad filter settings)
% output LaminarData contains the spike-phase distributions, circular
% resultants, and mean phase angles of the spike-LFP relationship on each
% channel.
% Example data
% exData1.mat = Marmoset Recording, exData2.mat = Macaque Recording, exData3.mat = Macaque Recording
% Load and run LaminarData = GetLaminarPhase(spikeTimes,LFP,chanMap,1000,[5 50],0);
%  10/21/2022 - ZD
%  02/23/2024 - HK


NChan = min(size(LFP));

% Flip the map so the top index is the most proximal electrode and the last
% index is the most distal electrode (i.e. deepest)
chanMap = flipud(chanMap);

% Initalize Vars
filtLFP = [];
phaseLFP = [];
phaseSpikes = cell(NChan,NChan);

% Filter and calculate LFP phase
% Chunk data because GP doesn't like big numbers
for i = 1:NChan
    for j = 1:NChan %For each spike channel
        if ~isempty(spikeTimes{chanMap(j)}) %No spikes? No phase.
            theseSpikes = spikeTimes{chanMap(j)};
            % Get rid of the spikes not in this segment
            theseSpikes = theseSpikes - min(thisSeg)+1;
            theseSpikes(theseSpikes < 1) = [];
            theseSpikes(theseSpikes > (Fs*10)) = [];
            %Get the spike phases
            phaseSpikes{i,j} = [phaseSpikes{i,j} phaseLFP(i,theseSpikes)];
        else
            phaseSpikes{i,j} = NaN;
        end
    end
end


spikePhase = [];
prefPhase = [];
for ii = 1:NChan
    for jj = 1:NChan
        spikePhase(ii,jj) = circ_r(phaseSpikes{jj,ii}'); %Calculate circular resultant
        prefPhase(ii,jj) = circ_mean(phaseSpikes{jj,ii}'); %Calculate mean phase angle
    end
end

phaseCorr = [];
for i = 1:NChan
    for j = 1:NChan
        phaseCorr(i,j) = circ_corrcc(phaseLFP(i,:),phaseLFP(j,:));
    end
end

%Collect Output
LaminarData.SPI = spikePhase;
LaminarData.Angle = prefPhase;
LaminarData.Dist = phaseSpikes;
LaminarData.PhaseCorr = phaseCorr;


figure
imagesc(spikePhase)
colormap hot
c = colorbar;
c.Label.String = 'SPI';
ylabel('Spiking Electrode Channel')
xlabel('LFP Electrode Channel')
set(gca,'fontsize',14,'linewidth',1.5)

figure
imagesc(prefPhase)
map = colorcet( 'C2' );
map = circshift(map,1);
colormap(map)
c = colorbar;
c.Label.String = 'Best Phase (rad)';
ylabel('Spiking Electrode Channel')
xlabel('LFP Electrode Channel')
set(gca,'fontsize',14,'linewidth',1.5)

figure
imagesc(phaseCorr)
colormap hot
c = colorbar;
c.Label.String = 'Circ Corr (r)';
ylabel('Channel')
xlabel('Channel')
set(gca,'fontsize',14,'linewidth',1.5)