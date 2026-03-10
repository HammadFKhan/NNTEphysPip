load('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\CCA_cooling\Day3DualShankCoolingRecording2.mat')
%%
cooledTrials = IntanBehaviour.hitTemp<-9;
baselineTrials = IntanBehaviour.hitTemp>=-9;
figure
subplot(1,4,[1 3])
Show_Spikes(M1Spikes.PSTH.hit.spks{12}),ylim([0 160]),xlim([1000,3000])  %1, 12, 13
subplot(1,4,4),plot(smoothdata(IntanBehaviour.hitTemp,1,'movmean',5),1:212),axis tight,ylim([0 160]) 
%%
cooledSpikes =binSpikes(M1Spikes.PSTH.hit.spks{12}(cooledTrials,:),20,20,1);
baseSpikes = binSpikes(M1Spikes.PSTH.hit.spks{12}(baselineTrials,:),20,20,1);
figure,hold on
plot(sum(baseSpikes,1)),plot(sum(cooledSpikes,1))
%%

figure
subplot(1,4,[1 3])
Show_Spikes(M2Spikes.PSTH.hit.spks{5}),ylim([0 160]),xlim([1000,3000]) %5
subplot(1,4,4),plot(smoothdata(IntanBehaviour.hitTemp,1,'movmean',5),1:212),axis tight,ylim([0 160]) 
%%
cooledSpikes =binSpikes(M2Spikes.PSTH.hit.spks{5}(cooledTrials,:),20,20,1);
baseSpikes = binSpikes(M2Spikes.PSTH.hit.spks{5}(baselineTrials,:),20,20,1);
figure,hold on
plot(mean(baseSpikes,1)),plot(mean(cooledSpikes,1))
%% Color code spikes
colorList=slanCM(100,150);
cd = [uint8((colorList)*255) uint8(ones(150,1))].';
n = 150;
%col = smoothdata(t,'movmean',10);
% Interp to make the line smoother
xin = interp1(1:150,x,1:0.05:150);
yin = interp1(1:150,y,1:0.05:150);
col = interp1(1:150,col,1:0.05:150);
col_map = (col - min(col))/(max(col)-min(col)) * (n-1) + 1;
for n = 1:length(col)
    cd1(:,n) = cd(:,floor(col_map(n)));
end
figure,
for n = 2:length(col)
    plot(xin(n-1:n),yin(n-1:n),'color',double(cd1(1:3,n))/255, 'LineWidth',2);hold on %cline( time, xf, [], angle(xgp) );
end
box off, axis off
%%
M = load( 'myMap.mat' );
smooth_temp = smoothdata(IntanBehaviour.hitTemp, 1, 'movmean', 10);
 
cmap = flip(M.myMap);
nColors = size(cmap,1);
y = 1:numel(smooth_temp);   % now vertical
v = smooth_temp(:);         % now horizontal
vMin = min(v);
vMax = max(v);
if vMax == vMin
    idx = ones(size(v));
else
    vNorm = (v - vMin) / (vMax - vMin);           % normalize to [0 1]
    idx = round( 1 + vNorm * (nColors - 1) );     % map to [1 nColors]
end
figure;
hold on;
for k = 1:(numel(y)-1)
    c = cmap(idx(k), :);
    % note the swapped order: x = smooth_temp, y = index
    plot(smooth_temp(k:k+1), y(k:k+1), 'Color', c, 'LineWidth', 1.5);
end
ylim([y(1) y(end)]);
xlabel('Temperature');
ylabel('Sample index');
title('Colored temperature trace with axes flipped');
colormap(cmap);
caxis([vMin vMax]);
colorbar;
ylim([0 160])

%%
function IntanBehaviour = grabTemp(IntanBehaviour,fpath)
if ~isfield(IntanBehaviour,'temperature')
    disp('No temp file added... correcting...')
    [filepath,~,~] = fileparts(fpath);
    load([filepath, '\loadme.mat']);
    if exist('ds_filename','var')
        data = matfile(ds_filename); % ds_filename comes from loadme.mat
    else
        data = matfile(ds_filename1); % ds_filename comes from loadme.mat
    end
    % check if data directory matches where the file originated; if not we note
    % the new directory path
    parameters.experiment = 'cue'; % self - internally generated, cue - cue initiated
    parameters.opto = 0; % 1 - opto ON , 0 - opto OFF
    parameters.cool = 1; % No Cool
    parameters.windowBeforePull = 1.5; % in seconds
    parameters.windowAfterPull = 1.5; % in seconds
    parameters.windowBeforeCue = 1.5; % in seconds
    parameters.windowAfterCue = 1.5; % in seconds
    parameters.windowBeforeMI = 1.5; % in seconds
    parameters.windowAfterMI = 1.5; % in seconds
    parameters.Fs = 1000; % Eventual downsampled data
    parameters.ts = 1/parameters.Fs;
    parameters.IntanFs = data.targetedFs;
    parameters.rows = 64;
    parameters.cols = 1;
    temperature = data.analogChannels(1,:);
    temperature = (temperature-1.25)/0.005;
    IntanBehaviour.temperature = resample(temperature,parameters.Fs,data.targetedFs);
    clear temperature
end
for n = 1:IntanBehaviour.nCueHit
    IntanBehaviour.hitTemp(n,1) = IntanBehaviour.temperature(IntanBehaviour.cueHitTrace(n).LFPIndex(1));
end
IntanBehaviour.hitTemp = IntanBehaviour.hitTemp-IntanBehaviour.temperature(100);
for n = 1:IntanBehaviour.nCueMiss
    IntanBehaviour.missTemp(n,1) = IntanBehaviour.temperature(IntanBehaviour.cueMissTrace(n).LFPIndex(1));
end
IntanBehaviour.missTemp = IntanBehaviour.missTemp-IntanBehaviour.temperature(100);
for n = 1:length(IntanBehaviour.missTrace)
    IntanBehaviour.FATemp(n,1) = IntanBehaviour.temperature(IntanBehaviour.missTrace(n).LFPIndex(1));
end
IntanBehaviour.FATemp = IntanBehaviour.FATemp-IntanBehaviour.temperature(100);
end

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
