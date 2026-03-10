%% Code to prep the coupling data for TW and spikes
figure,
t = mean(waveDynamics.rawWaveSpeedhit)';
t = t(1:20:end-1);
col = [smoothdata((t),'gaussian',10)'];
plotNeuralTrajWave(neuralDynamics.hit.r(1,:),neuralDynamics.hit.r(2,:),col)
%% Do stats on the trajectory and wave coupling
close all
[wavePGDCoupling,waveSpeedCoupling] = getTrajectoryWaveStats(neuralDynamics,waveDynamics);
%%
f = figure;
f.Position = [1607 467 560 820];
subplot(231),plot(squeeze(M1neuralDynamics.hit.X(1,:,:)),'color',[0.5 0.5 0.5 0.1]);hold on,plot(mean(squeeze(M1neuralDynamics.hit.X(1,:,:)),2),'color',[1 0.5 0.5]);
subplot(234),plot(squeeze(M1neuralDynamics.hit.X(2,:,100:end)),'color',[0.5 0.5 0.5 0.1]);hold on,plot(mean(squeeze(M1neuralDynamics.hit.X(2,:,100:end)),2),'color',[1 0.5 0.5]);

subplot(232),plot(squeeze(M1neuralDynamics.miss.X(1,:,:)),'color',[0.5 0.5 0.5 0.1]);hold on,plot(mean(squeeze(M1neuralDynamics.miss.X(1,:,:)),2),'color',[1 0.5 0.5]);
subplot(235),plot(squeeze(M1neuralDynamics.miss.X(2,:,1:end)),'color',[0.5 0.5 0.5 0.1]);hold on,plot(mean(squeeze(M1neuralDynamics.miss.X(2,:,1:end)),2),'color',[1 0.5 0.5])

subplot(233),plot(squeeze(M1neuralDynamics.MIFA.X(1,:,:)),'color',[0.5 0.5 0.5 0.1]);hold on,plot(mean(squeeze(M1neuralDynamics.MIFA.X(1,:,:)),2),'color',[1 0.5 0.5]);
subplot(236),plot(squeeze(M1neuralDynamics.MIFA.X(2,:,1:end)),'color',[0.5 0.5 0.5 0.1]);hold on,plot(mean(squeeze(M1neuralDynamics.MIFA.X(2,:,1:end)),2),'color',[1 0.5 0.5])
%%
close all
for neuron = 1:20;
f = figure(neuron);
clf
f.Position = [1607 467 560 820];
subplot(3,1,1)
Show_Spikes(Spikes.PSTH.hit.spks{neuron}),box off,set(gca,'xtick',[],'tickdir','out') ;
subplot(3,1,2)
Show_Spikes(Spikes.PSTH.miss.spks{neuron}),box off,set(gca,'xtick',[],'tickdir','out');
subplot(3,1,3),plot(mean(binSpikes(Spikes.PSTH.hit.spks{neuron},20,20,1))*20);hold on
plot(mean(binSpikes(Spikes.PSTH.miss.spks{neuron},20,20,1))*20);
end

for neuron = 1:20;
f = figure(20+neuron);
clf
f.Position = [1607 467 560 820];
subplot(3,1,1)
Show_Spikes(Spikes.PSTH.MIHit.spks{neuron}),box off,set(gca,'xtick',[],'tickdir','out') ;
subplot(3,1,2)
Show_Spikes(Spikes.PSTH.MIFA.spks{neuron}),box off,set(gca,'xtick',[],'tickdir','out');
subplot(3,1,3),plot(mean(binSpikes(Spikes.PSTH.MIHit.spks{neuron},20,20,1))*20);hold on
plot(mean(binSpikes(Spikes.PSTH.MIFA.spks{neuron},20,20,1))*20);
end
%%
PGD = M1Waves
%%
function [binned_smoothed, t_bins] = binSpikes(spikes, bin_size, smooth_sigma, dt)
% BIN_AND_SMOOTH_SPIKES
%   Bin and Gaussian-smooth spike data (trials x timebins).
%
%   Inputs
%   ------
%   spikes       : [nTrials x nTime] matrix of spikes (0/1 or counts)
%   bin_size     : integer, number of original samples per bin (>= 1)
%   smooth_sigma : std of Gaussian kernel in same time units as dt
%                  (if 0 or empty, no smoothing)
%   dt           : sampling interval of 'spikes' (time units per sample)
%
%   Outputs
%   -------
%   binned_smoothed : [nTrials x nBins] binned and smoothed counts
%   t_bins          : [1 x nBins] bin-center times (same units as dt)

    if nargin < 4 || isempty(dt)
        dt = 1;
    end
    if nargin < 3 || isempty(smooth_sigma)
        smooth_sigma = 0;
    end
    if nargin < 2 || isempty(bin_size)
        bin_size = 1;
    end

    if ndims(spikes) ~= 2
        error('spikes must be a 2D (trials x time) array.');
    end

    [nTrials, nTime] = size(spikes);

    if bin_size < 1
        error('bin_size must be >= 1.');
    end
    if bin_size > nTime
        error('bin_size cannot exceed number of timepoints.');
    end

    % ----- Bin along time axis -----
    nBins   = floor(nTime / bin_size);
    trimmed = spikes(:, 1:nBins * bin_size);

    % reshape: [nTrials x nBins x bin_size], then sum over bin_size
    trimmed_reshaped = reshape(trimmed, nTrials, bin_size, nBins);
    binned = squeeze(sum(trimmed_reshaped, 2));   % [nTrials x nBins]

    % Time axis for bin centers
    bin_dt = bin_size * dt;
    t_bins = ((0:nBins-1) + 0.5) * bin_dt;

    % ----- Gaussian smoothing along time axis -----
    if smooth_sigma > 0
        % sigma in bins
        sigma_bins = smooth_sigma / bin_dt;

        % Kernel extent: +/- 3 sigma
        half_width = ceil(3 * sigma_bins);
        x = -half_width:half_width;

        % Gaussian kernel
        gauss = exp(-0.5 * (x ./ sigma_bins).^2);
        gauss = gauss / sum(gauss);  % normalize to area 1

        % Convolve each trial with same kernel along time (conv2)
        % conv2 size option 'same' keeps [nTrials x nBins]
        binned_smoothed = conv2(binned, gauss, 'same');
    else
        binned_smoothed = binned;
    end
end