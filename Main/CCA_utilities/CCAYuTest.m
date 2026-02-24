%% Plot out spike data
[~] = sortSpkLever(M1Spikes,IntanBehaviour);
[~] = sortSpkLever(M2Spikes,IntanBehaviour);
%%
M1SpikesCCA = [];
M2SpikesCCA = [];
for n = 1:length(M1Spikes.GPFA.hit.dat)
    M1SpikesCCA(:,:,n) = M1Spikes.GPFA.hit.dat(n).spikes;
    M2SpikesCCA(:,:,n) = M2Spikes.GPFA.hit.dat(n).spikes;
end
figure,subplot(121),imagesc(sum(M1SpikesCCA,3))
subplot(122),imagesc(sum(M2SpikesCCA,3))

%% CCA Boosting
% M1, M2: p1 x T x N and p2 x T x N, binary 0/1 spike trains at 1 ms
[p1, T, N] = size(M1SpikesCCA);
[p2, ~, ~] = size(M2SpikesCCA);
winStart = 1500;         % start time bin (inclusive, 1-based)
winEnd   = 1780;         % end time bin (inclusive)
nSharedNeur1 = min(p1,15);       % number of neurons in M1 to carry shared pattern
nSharedNeur2 = min(p2,15);       % number of neurons in M2 to carry shared pattern
trialFrac    = 1.0;     % fraction of trials to modify (0–1)
rateBoost    = 0.001;     % probability of setting a spike at shared times

% choose which trials to boost
nBoostTrials    = round(trialFrac * N);
boostTrialIndex = randperm(N, nBoostTrials);

% copy originals
M1_mod = M1SpikesCCA;
M2_mod = M2SpikesCCA;

% choose neurons that will carry shared pattern
neurIdx1 = randperm(p1, nSharedNeur1);
neurIdx2 = randperm(p2, nSharedNeur2);

for n = boostTrialIndex
    % define a shared binary pattern over time in the window
    % (same across all chosen neurons in both areas, for this trial)
    shared_pattern = rand(1, winEnd-winStart+1) < rateBoost;  % 0/1

    % add this pattern to selected neurons in M1 and M2
    for i = 1:nSharedNeur1
        m1n = neurIdx1(i);
        % OR with existing spikes to keep binary
        M1_mod(m1n, winStart:winEnd, n) = ...
            M1_mod(m1n, winStart:winEnd, n) | shared_pattern;
    end

    for j = 1:nSharedNeur2
        m2n = neurIdx2(j);
        M2_mod(m2n, winStart:winEnd, n) = ...
            M2_mod(m2n, winStart:winEnd, n) | shared_pattern;
    end
end

figure,subplot(121),imagesc(sum(M1_mod,3))
subplot(122),imagesc(sum(M2_mod,3))
M1SpikesCCA = M1_mod;
M2SpikesCCA = M2_mod;
%% save new data
for n = 1:size(M1SpikesCCA,3)
    M1Spikes.GPFA.hit.dat(n).spikes = M1SpikesCCA(:,:,n);
    M2Spikes.GPFA.hit.dat(n).spikes = M2SpikesCCA(:,:,n);
end
spath = 'Y:\Hammad\Ephys\LeverTask\Data_for_Figures\M1M2DualShank\CCA';
[path,name,ext] = fileparts(fpath);
sessionName = [spath,'\',name(1:end-14)];

save(sessionName,"M2Spikes","M1Spikes","IntanBehaviour","fpath","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
disp('Saved CCA data!')
%% Byron's code
% To avoid including a large data file in the repository, the data
% loaded here have been pre-binned using 10ms windows, i.e., each entry in
% the spikes matrices indicates the number of spikes recorded in a 10ms
% window
% Wrap spike arrays into Byron's format
spikes = cell(2,1);
spikes{1} = M1SpikesCCA;   % p1 x T x N
spikes{2} = M2SpikesCCA;   % p2 x T x N

N = size(M1SpikesCCA,3);
expCond = ones(N,1);   % or zeros(N,1); label 1 for all trials

%% Example computation of a cross-correlation map
% Figs. 3 and 4

% The units of the arguments are with respect to the binning window used
% to bin spikes.

argIn.BinWidth = 10;     % 10ms
argIn.MaxDelay = 20;    % 100ms
argIn.TimeStep = 40;     % 40ms
argIn.WindowLength = 80; % 80ms

%argIn.NumWorkers = Inf; % Requires Parallel Processing Toolbox

disp(argIn)
argOut = ComputeCorrMap(spikes, expCond, argIn);

%%
CANONICAL_PAIR_IDX = 1;
mapDim = size(argOut.CorrMap, 2);
delays = (-argIn.MaxDelay:argIn.MaxDelay)*10; % Convert to ms
t = (1:argIn.TimeStep:argIn.TimeStep*mapDim)*10; % Convert to ms

figure(1);

imagesc( delays, t, argOut.CorrMap(:,:,CANONICAL_PAIR_IDX)' )

ax = gca;
ax.YDir = 'Normal';

xlabel('Delay')
ylabel('Time')

figure,plot(t,smoothdata(squeeze(argOut.CorrMap(21,:,CANONICAL_PAIR_IDX))),'linewidth',2)
xlabel('Time')
%% Example computation of the interaction structure analysis
% Fig. 6

clear argIn

argIn.TimePeriods = [...
    (  0:20:40)' ( 20:20:60)'; ...
    (128:20:168)' (148:20:188)'] + 5;

% Can take up to 15min due to the 10-fold cross-validation
argOut = CovStabilityAcrossTimeAnalysis(spikes, expCond, argIn);

%%
RANK_TO_PLOT = 2;

normFactor = diag(argOut.CvR(:,:,RANK_TO_PLOT));
numTimePeriods = size(argIn.TimePeriods, 1);

figure(2);

imagesc(argOut.CvR(:,:,RANK_TO_PLOT)./repmat(normFactor', numTimePeriods, 1))

ax = gca;
ax.YDir = 'Normal';

axis square

xlabel('Time Used For Correlation')
ylabel('Time Used For Fitting')


