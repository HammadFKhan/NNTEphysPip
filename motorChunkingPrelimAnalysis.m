%%
files = dir(fullfile('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\M1_GSP','*.mat'));
M1neuralDynamics = struct();
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    if ~isfield(IntanBehaviour, 'parameters')
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
        parameters.IntanFs = 2000;
        parameters.rows = 64;
        parameters.cols = 1;
        IntanBehaviour.parameters = parameters;
    end
    Spikes = makeSpikeGPFA(Spikes);
    Spikes.GPFA.HitMiss.dat = [Spikes.GPFA.hit.dat,Spikes.GPFA.miss.dat];
    for n = 1:IntanBehaviour.nCueHit%+1:length(Spikes.GPFA.BaselineOpto.dat) %fix trials
        Spikes.GPFA.HitMiss.dat(n).trialId = n;
    end
    Spikes.GPFA.MIHitFA.dat = [Spikes.GPFA.MIHit.dat,Spikes.GPFA.MIFA.dat];
    for n = length(IntanBehaviour.MIHitTrace)+1:length(Spikes.GPFA.MIHitFA.dat) %fix trials
        Spikes.GPFA.MIHitFA.dat(n).trialId = n;
    end
    addpath(genpath('C:\Users\khan332\Documents\GitHub\Ephy2\NeuralTraj'));
    addpath(genpath('mat_results'));
    if exist('mat_results','dir'),rmdir('mat_results','s'),end
    [Spikes.GPFA.resultHit,Spikes.GPFA.seqTrainHit] = gpfaAnalysis(Spikes.GPFA.hit.dat,1); %Run index
    [Spikes.GPFA.resultMiss,Spikes.GPFA.seqTrainMiss] = gpfaAnalysis(Spikes.GPFA.miss.dat,2); %Run index
    [Spikes.GPFA.resultMIHit,Spikes.GPFA.seqTrainMIHit] = gpfaAnalysis(Spikes.GPFA.MIHit.dat,3); %Run index
    [Spikes.GPFA.resultMIFA,Spikes.GPFA.seqTrainMIFA] = gpfaAnalysis(Spikes.GPFA.MIFA.dat,4); %Run index
    [Spikes.GPFA.resultHitMiss,Spikes.GPFA.seqTrainHitMiss] = gpfaAnalysis(Spikes.GPFA.HitMiss.dat,5); %Run index
    [Spikes.GPFA.resultMIHitFA,Spikes.GPFA.seqTrainMIHitFA] = gpfaAnalysis(Spikes.GPFA.MIHitFA.dat,6); %Run index
    M1neuralDynamics(fileNum).fname = files(fileNum).name;
    M1neuralDynamics(fileNum).IntanBehaviour = IntanBehaviour;
    [M1neuralDynamics(fileNum).neuralDynamics,M1waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);
    M1neuralDynamics(fileNum).Spikes = Spikes;
    close all
end
%% Preliminary analysis on motor chunking in mice
leverTot = [];
for n = 1:length(M1neuralDynamics)
    IntanBehaviour = M1neuralDynamics(n).IntanBehaviour ;
    %%%
    leverTrace = arrayfun(@(x) x.trace',IntanBehaviour.cueHitTrace,'uniformoutput',0);
    leverTrace = vertcat(leverTrace{:});
    leverTrace = (leverTrace-min(leverTrace,[],2))./(max(leverTrace,[],2)-min(leverTrace,[],2));
    leverTot = vertcat(leverTot,leverTrace);
end
leverTrace = leverTot;
reconstructedData = getPCALever(leverTrace);
leverTrace = reconstructedData;
%%
num_trials = size(leverTrace, 1);
double_pull_trials = [];

peaksTot = [];
locsTot = [];
for trial = 1:num_trials
    [peaks, locs] = findpeaks(leverTrace(trial, 1500:end), 'MinPeakProminence', 0.50, 'MinPeakDistance', 100);
    if length(peaks) >= 2
        double_pull_trials = [double_pull_trials, trial];
        peaksTot(trial,:) = peaks(1:2);
        locsTot(trial,:) = locs(1:2)+1500;
    end
end
id = find(peaksTot(:,1)==0);
peaksTot(id,:) = [];
locsTot(id,:) = [];

figure(3),clf;
for i = 1:9
    trial = double_pull_trials(i);
    subplot(3,3, i);
    plot(leverTrace(trial, :));
    hold on;
    [peaks, locs] = findpeaks(leverTrace(trial, 1500:end), 'MinPeakProminence', 0.50, 'MinPeakDistance', 100);
    plot(locs+1500, peaks, 'ro');
    title(['Trial ', num2str(trial)]);
    xlabel('Time');
    ylabel('Lever Position');
end
%% show example single vs compound traces
figure,
subplot(211),plot(smoothdata(leverTrace(8,:),'movmean',50),'k','linewidth',2)
box off, axis square
subplot(212),plot(smoothdata(leverTrace(4,:),'movmean',50),'linewidth',2,'color',[0.5 0.5 0.5])
box off, axis square
%%
pullITI = locsTot(:,2)-locsTot(:,1);
figure,
subplot(121),histogram(pullITI,0:100:1500,'normalization','probability','edgecolor','none')
box off, axis square
set(gca,'tickdir','out','fontsize',12),xlabel('Time between pulls (ms)'),ylabel('Probability')
subplot(122),histogram(peaksTot(:,1)-peaksTot(:,2),-0.7:0.1:0.7,'normalization','probability','edgecolor','none')
box off, axis square
set(gca,'tickdir','out','fontsize',12),xlabel('Amplitude difference from first pull'),ylabel('Probability')
%% show combined single vs compound traces
motorChunkshort = find(pullITI<600);
singlePullTrials = 1:size(leverTrace,1);
singlePullTrials(double_pull_trials) = [];
figure,
subplot(311),plot(mean(smoothdata(leverTrace(singlePullTrials,:),2,'movmean',50)),'k','linewidth',3)
box off, axis square,xlim([1000 3000])
subplot(312),plot(mean(smoothdata(leverTrace(double_pull_trials(motorChunkshort),:),2,'movmean',50)),'linewidth',3,'color',[0.3 0.3 0.3])
box off, axis square,xlim([1000 3000])
motorChunklong = find(pullITI>1000);
subplot(313),plot(mean(smoothdata(leverTrace(double_pull_trials(motorChunklong),:),2,'movmean',50)),'linewidth',3,'color',[0.6 0.6 0.6])
box off, axis square,xlim([1000 3000])


%% Build neural dynamics across all mice to find compound lever movement

% Append data to the trials 
latentDynamics = [];
for n = 1:length(M1neuralDynamics)
    latentDynamics = cat(3,latentDynamics,M1neuralDynamics(n).neuralDynamics.hit.X);
end

% Seperate latent dynamics based on behaviour trials
% assert that latent dynamic trials is equal to lever trials
assert(size(latentDynamics,3)==size(leverTrace,1))
latentDynamicsSingle = latentDynamics(:,:,singlePullTrials);
latentDynamicsChunkShort = latentDynamics(:,:,double_pull_trials(motorChunkshort));
latentDynamicsChunkLong = latentDynamics(:,:,double_pull_trials(motorChunklong));

%%% Plot it out
figure,
subplot(311),plot(mean(squeeze(latentDynamicsSingle(1,:,:)),2),'linewidth',2,'color',[0 0 0])
subplot(312),plot(mean(squeeze(latentDynamicsSingle(2,:,:)),2),'linewidth',2,'color',[0.3 0.3 0.3])
subplot(313),plot(mean(squeeze(latentDynamicsSingle(3,:,:)),2),'linewidth',2,'color',[0.6 0.6 0.6])

figure,
subplot(311),plot(mean(squeeze(latentDynamicsChunkShort(1,:,:)),2),'linewidth',2,'color',[0 0 0])
subplot(312),plot(mean(squeeze(latentDynamicsChunkShort(2,:,:)),2),'linewidth',2,'color',[0.3 0.3 0.3])
subplot(313),plot(mean(squeeze(latentDynamicsChunkShort(3,:,:)),2),'linewidth',2,'color',[0.6 0.6 0.6])

figure,
subplot(311),plot(mean(squeeze(latentDynamicsChunkLong(1,:,:)),2),'linewidth',2,'color',[0 0 0])
subplot(312),plot(mean(squeeze(latentDynamicsChunkLong(2,:,:)),2),'linewidth',2,'color',[0.3 0.3 0.3])
subplot(313),plot(mean(squeeze(latentDynamicsChunkLong(3,:,:)),2),'linewidth',2,'color',[0.6 0.6 0.6])

%% LOCAL FUNCTIONS
function reconstructedData = getPCALever(data)
% Transpose the data
data = data';

% Perform PCA
[coeff, score, latent] = pca(data);

% Determine how many principal components to retain
cum_var = cumsum(latent ./ sum(latent));
n_components = find(cum_var >= 0.95, 1, 'first');

% Select the first n_components principal components
selected_coeff = coeff(:, 1:n_components);
selected_score = score(:, 1:n_components);

% Calculate the mean of your original data
mean_data = mean(data);

% Reconstruct the data using the selected principal components
reconstructedData = selected_score * selected_coeff' + repmat(mean_data, size(data, 1), 1);
reconstructedData = reconstructedData';
end
