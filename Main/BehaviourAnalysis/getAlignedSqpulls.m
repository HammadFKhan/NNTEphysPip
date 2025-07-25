function warpedSpks = getAlignedSqpulls(Spikes,warpedSpks,IntanBehaviour)

%% Plot trial aligned by first pull
% than by earliest second pull

pullIndex = floor([warpedSpks.pull1A.pull1',warpedSpks.pull1A.pull2',warpedSpks.pull1A.pull3']);
sqPullMask1 = getLeverplots(pullIndex,IntanBehaviour);

% Plot by second pull
pullIndex = floor([warpedSpks.pull2A.pull1',warpedSpks.pull2A.pull2',warpedSpks.pull2A.pull3']);
sqPullMask2 = getLeverplots(pullIndex,IntanBehaviour);

% Plot by second pull
pullIndex = floor([warpedSpks.pull3A.pull1',warpedSpks.pull3A.pull2',warpedSpks.pull3A.pull3']);
sqPullMask3 = getLeverplots(pullIndex,IntanBehaviour);


%% plot the warp line shift
bin = 10;
color = [46,49,179;46,179,49;179,49,46]/255;
figure,
binned_data = getBin(sqPullMask1,bin);
plot(smoothdata(sum(binned_data),'gaussian',50),'color',color(1,:)),hold on
binned_data = getBin(sqPullMask2,bin);
plot(smoothdata(sum(binned_data),'gaussian',50),'color',color(2,:)),hold on
binned_data = getBin(sqPullMask3,bin);
plot(smoothdata(sum(binned_data),'gaussian',50),'color',color(3,:)),hold on
box off,set(gca,'tickdir','out'),axis square
xlim([100 400])

%% sort to first pull
% Nice neurons 1 3 6 7 10 19 17


% [warpSpikes,warpTime] = getwarpedSpikes(warpedSpks.pull3A.Spks,warpedSpks.pull3A.pull1,IntanBehaviour);
% figure
% count = 1;
% color = [46,49 149]/255;
% for neuron =  1:30
%     subplot(6,5,count)
%     ylim([1 100])
%     xline(0,'k','reward')
%     spkTemp = squeeze(warpSpikes(:,:,neuron));
%     spkTemp(spkTemp==0) = NaN;
%     for n = 1:size(spkTemp,1)
%         scatter(warpTime,n*spkTemp(n,:),5,'filled','MarkerFaceColor',[0,0,0]),hold on
%     end
%     ylim([1 100])
%     xline(0,'k','reward')
%     xlim([warpTime(1), warpTime(end)])
%     title(['Neuron ' num2str(neuron)])
%     count = count+1;
% %     for n = 1:size(trialMask,1)
% %         scatter(time(ft(n)),n*trialMask(firstPull(n),ft(n)),10,'filled','MarkerEdgeColor',color,'MarkerFaceColor','none'),hold on
% %     end
% end
%% Plot PSTH of example neurons with shifted responses
% Neurons to look at [3 4 23 24]
f = figure;tiledlayout(2,4)
f.Position = [680 558 960 320];

count = 1;
[warpSpikes,warpTime] = getwarpedSpikes(warpedSpks.pull3A.Spks,warpedSpks.pull3A.pull1,IntanBehaviour);
for neuron = [3 4 23 24]
    spkTemp = squeeze(warpSpikes(:,:,neuron));
    nexttile
    spkTemp(spkTemp==0) = NaN;
    for n = 1:size(spkTemp,1)
        scatter(warpTime,n*spkTemp(n,:),0.5,'filled','MarkerFaceColor',[0,0,0]),hold on
    end
    ylim([0 100])
    %xline(0,'k')
    xlim([warpTime(1), warpTime(end)])
    xlabel('Time from pull (s)')
    title(['Neuron ' num2str(neuron)])
end
for neuron = [3 4 23 24]
    spkTemp = squeeze(warpSpikes(:,:,neuron));
    nexttile
    bin = 20;
    binnedSpk= getBin(spkTemp,bin);
    plot(linspace(warpTime(1), warpTime(end),length(binnedSpk)),smoothdata(sum(binnedSpk)*(1000/bin),'gaussian',20))
    xlim([warpTime(1), warpTime(end)])
    xlabel('Time from pull (s)'),box off,set(gca,'tickdir','out')
end
%%
f = figure;tiledlayout(1,4)
f.Position = [680 558 960 320];
color = [46,49,179,55;46,179,49,55;179,49,46,55]/255;
nexttile
[warpSpikes,warpTime] = getwarpedSpikes(warpedSpks.pull1A.Spks,warpedSpks.pull1A.pull1,IntanBehaviour);
warpedSpks.pull1A.warpSpikes = warpSpikes;
warpedSpks.pull1A.warpTime = warpTime;
binnedSpk = [];

bin = 20;
for neuron = 1:Spikes.nSpikes
spkTemp = squeeze(warpSpikes(:,:,neuron));
binnedSpk(neuron,:) = zscore(sum(getBin(spkTemp,bin))*(1000/bin),1);
end
for n = 1:Spikes.nSpikes
    plot(linspace(warpTime(1), warpTime(end),length(binnedSpk)),smoothdata(binnedSpk(n,:),'gaussian',20),'color',color(1,:)),hold on %smoothdata(zscore(binnedSpk,0,2)','gaussian',1))
end
    plot(linspace(warpTime(1), warpTime(end),length(binnedSpk)),smoothdata(mean(binnedSpk),'gaussian',20),'k'),hold on %smoothdata(zscore(binnedSpk,0,2)','gaussian',1))

xlim([warpTime(1), warpTime(end)])
xlabel('Time from pull (s)'),box off,set(gca,'tickdir','out'),axis square
ylim([-2.5 3.5])

nexttile
[warpSpikes,warpTime] = getwarpedSpikes(warpedSpks.pull2A.Spks,warpedSpks.pull2A.pull2,IntanBehaviour);
warpedSpks.pull2A.warpSpikes = warpSpikes;
warpedSpks.pull2A.warpTime = warpTime;

bin = 20;
for neuron = 1:Spikes.nSpikes
spkTemp = squeeze(warpSpikes(:,:,neuron));
binnedSpk(neuron,:) = zscore(sum(getBin(spkTemp,bin))*(1000/bin),1);
end
for n = 1:Spikes.nSpikes
    plot(linspace(warpTime(1), warpTime(end),length(binnedSpk)),smoothdata(binnedSpk(n,:),'gaussian',20),'color',color(2,:)),hold on %smoothdata(zscore(binnedSpk,0,2)','gaussian',1))
end
    plot(linspace(warpTime(1), warpTime(end),length(binnedSpk)),smoothdata(mean(binnedSpk),'gaussian',20),'k'),hold on %smoothdata(zscore(binnedSpk,0,2)','gaussian',1))

xlim([warpTime(1), warpTime(end)])
xlabel('Time from pull (s)'),box off,set(gca,'tickdir','out'),axis square
ylim([-2.5 3.5])

nexttile
[warpSpikes,warpTime] = getwarpedSpikes(warpedSpks.pull3A.Spks,warpedSpks.pull3A.pull3,IntanBehaviour);
warpedSpks.pull3A.warpSpikes = warpSpikes;
warpedSpks.pull3A.warpTime = warpTime;

bin = 20;
for neuron = 1:Spikes.nSpikes
spkTemp = squeeze(warpSpikes(:,:,neuron));
binnedSpk(neuron,:) = zscore(sum(getBin(spkTemp,bin))*(1000/bin),1);
end
for n = 1:Spikes.nSpikes
    plot(linspace(warpTime(1), warpTime(end),length(binnedSpk)),smoothdata(binnedSpk(n,:),'gaussian',20),'color',color(3,:)),hold on %smoothdata(zscore(binnedSpk,0,2)','gaussian',1))
end
    plot(linspace(warpTime(1), warpTime(end),length(binnedSpk)),smoothdata(mean(binnedSpk),'gaussian',20),'k'),hold on %smoothdata(zscore(binnedSpk,0,2)','gaussian',1))

xlim([warpTime(1), warpTime(end)])
xlabel('Time from pull (s)'),box off,set(gca,'tickdir','out'),axis square
ylim([-2.5 3.5])

nexttile
[warpSpikes,warpTime] = getwarpedSpikes(warpedSpks.pull23A.Spks,warpedSpks.pull23A.pull1,IntanBehaviour);
bin = 20;
for neuron = 1:Spikes.nSpikes
spkTemp = squeeze(warpSpikes(:,:,neuron));
binnedSpk(neuron,:) = zscore(sum(getBin(spkTemp,bin))*(1000/bin),1);
end
for n = 1:Spikes.nSpikes
    plot(linspace(warpTime(1), warpTime(end),length(binnedSpk)),smoothdata(binnedSpk(n,:),'gaussian',20),'color',[0.5 0.5 0.5 0.25]),hold on %smoothdata(zscore(binnedSpk,0,2)','gaussian',1))
end
    plot(linspace(warpTime(1), warpTime(end),length(binnedSpk)),smoothdata(mean(binnedSpk),'gaussian',20),'k'),hold on %smoothdata(zscore(binnedSpk,0,2)','gaussian',1))

xlim([warpTime(1), warpTime(end)])
xlabel('Time from pull (s)'),box off,set(gca,'tickdir','out'),axis square
ylim([-2.5 3.5])

%% OUTPUT full warpedSpks and times
warpSpikestot = zeros([size(warpSpikes),3]);
warpSpikestot(:,:,:,1) = warpedSpks.pull1A.warpSpikes;
warpSpikestot(:,:,:,2) = warpedSpks.pull2A.warpSpikes;
warpSpikestot(:,:,:,3) = warpedSpks.pull3A.warpSpikes;
warpedSpks.warpSpikes = warpSpikestot;

warpedTimetot(:,1) = warpedSpks.pull1A.warpTime;
warpedTimetot(:,2) = warpedSpks.pull2A.warpTime;
warpedTimetot(:,3) = warpedSpks.pull3A.warpTime;
warpedSpks.warpedTime = warpedTimetot;

end
%% LOCAL FUNCTIONS
function sqPullMask = getLeverplots(pullIndex,IntanBehaviour)

sqPullMask = zeros(length(IntanBehaviour.hitTrace),length(IntanBehaviour.hitTrace(1).time));
for n = 1:length(IntanBehaviour.hitTrace)
    sqPullMask(n,pullIndex(n,:)) = 1;
end

% Plot out interleaved flagging of the contextual cue data
% Return interleaved trials index of when each cue was initiated
trialMask = sqPullMask;
trialMask(trialMask==0) = NaN;

warpTimeWin = [-IntanBehaviour.parameters.windowBeforePull,IntanBehaviour.parameters.windowAfterPull];
time = linspace(warpTimeWin(1),warpTimeWin(2),size(trialMask,2));
warpTime = time-time(median(pullIndex(:,1)));
figure
color = [46,49,149]/255;
for n = 1:size(trialMask,1)
    scatter(warpTime,n*trialMask(n,:),10,'filled','MarkerFaceColor',color),hold on
end
xlim([warpTime(1) warpTime(end)])
ylim([1 100])
axis square
end

function [warpSpikes,warpTime] = getwarpedSpikes(spkTemp,pullTemp,IntanBehaviour)

% Define time boundaries
t_min = 0;
t_max = length(IntanBehaviour.hitTrace(1).time);

% Logical index of spikes within bounds
valid_idx = (spkTemp.spiketimes >= t_min) & (spkTemp.spiketimes <= t_max);

% Keep only valid data
trial_ids = spkTemp.trials(valid_idx);
neuron_ids = spkTemp.neurons(valid_idx);
spiketimes = spkTemp.spiketimes(valid_idx);
% Convert to 1-based indexing
trial_ids_1based = trial_ids + 1;
neuron_ids_1based = neuron_ids + 1;
spiketimes_1based = floor(spiketimes) + 1;


% Dimensions (assuming trial_ids, spiketimes, neuron_ids all 1-based)
nTrials = max(trial_ids_1based);
nTimeBins = t_max;
nNeurons = max(neuron_ids_1based);

% Initialize the 3D matrix with dimension: trials x time x neurons
warpSpikes = zeros(nTrials, nTimeBins, nNeurons);

% Compute linear indices with the appropriate order of subscripts
idx = sub2ind([nTrials, nTimeBins, nNeurons], trial_ids_1based, spiketimes_1based, neuron_ids_1based);

% Assign spikes
warpSpikes(idx) = 1;

warpTimeWin = [-IntanBehaviour.parameters.windowBeforePull,IntanBehaviour.parameters.windowAfterPull];
time = linspace(warpTimeWin(1),warpTimeWin(2),t_max);
warpTime = time-time(floor(median(pullTemp)));
end

function binned_data = getBin(data,bin_size)
[n_trials, n_time] = size(data);

% Trim time dimension to multiple of bin_size
n_time_trim = floor(n_time / bin_size) * bin_size;
data_trim = data(:, 1:n_time_trim);

% Reshape to [n_trials, bin_size, n_time_trim/bin_size]
data_reshaped = reshape(data_trim', bin_size, [], n_trials); 
% Note the transpose is so time is first dimension for reshaping

% Sum or average within bins (along first dimension)
binned_data = squeeze(mean(data_reshaped, 1))';  
% Output size: [n_trials, n_time_trim/bin_size]
end