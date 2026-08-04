%% Bipoles in M2 recording in M1; time from first spike latency as a function of opto pulsing. 
% First reconstructu opto traces per hit trial
% After window we then extract the opto cue
files = dir(fullfile('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\M2BiPOLES\excitationM1R\','*.mat'));
firstSpikeBiPoleTotal = {};
taggedSpkTotal = {};
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    optoTraceHit = [];
    for n = 1:length(IntanBehaviour.optoCueHitTrace)
        win = [IntanBehaviour.optoCueHitTrace(n).LFPIndex(1), IntanBehaviour.optoCueHitTrace(n).LFPIndex(end)];
        optoTraceHit(n,:) = IntanBehaviour.optoPulseTrace(win(1):win(2));
    end
    % Now calculate time from first spike based on opto tag M1 neurons and the
    % time in which the pulse reaches one.
    % BiPOLES index is based off the hit PSTH
    optoTrials = vertcat(IntanBehaviour.cueHitTrace.opto);   % 0/1 labels
    optoTrials = find(optoTrials==1);
    % Do some assertions so we know we got it right
    assert(size(optoTrials,1)==size(optoTraceHit,1))
    assert(size(vertcat(IntanBehaviour.cueHitTrace.opto),1)==size(Spikes.PSTH.hit.spks{1},1))
    count = 1;
    firstSpikeBiPole = []
    for opto = 1:size(optoTraceHit,1)
        optoId = find(diff(optoTraceHit(opto,:))==1)+1;
        for neuron = find(Spikes.BiPOLES.tagged==1)
            disp(['Calcuating neuron ', num2str(neuron) '...'])
            BiPolesNeuron = Spikes.PSTH.hit.spks{neuron};
            for optoIdx = 1:length(optoId)
                if optoIdx==length(optoId)
                    firstSpk = find(BiPolesNeuron(optoTrials(opto),optoId(optoIdx):end)==1,1);
                else
                    firstSpk = find(BiPolesNeuron(optoTrials(opto),optoId(optoIdx):(optoId(optoIdx+1)-10))==1,1);
                end
                if ~isempty(firstSpk)
                    firstSpikeBiPole(count) = firstSpk;
                    count = count+1;
                end
            end
        end
    end
    firstSpikeBiPoleTotal{fileNum} = firstSpikeBiPole;
    taggedSpkTotal{fileNum} = Spikes.BiPOLES.tagged;

    % now grab the non tagged
    optoTraceHit = [];
    for n = 1:length(IntanBehaviour.cueHitTrace)
        win = [IntanBehaviour.optoCueHitTrace(n).LFPIndex(1), IntanBehaviour.optoCueHitTrace(n).LFPIndex(end)];
        optoTraceHit(n,:) = IntanBehaviour.optoPulseTrace(win(1):win(2));
    end
    % Now calculate time from first spike based on opto tag M1 neurons and the
    % time in which the pulse reaches one.
    % BiPOLES index is based off the hit PSTH
    optoTrials = vertcat(IntanBehaviour.cueHitTrace.opto);   % 0/1 labels
    optoTrials = find(optoTrials==1);
    % Do some assertions so we know we got it right
    assert(size(optoTrials,1)==size(optoTraceHit,1))
    assert(size(vertcat(IntanBehaviour.cueHitTrace.opto),1)==size(Spikes.PSTH.hit.spks{1},1))
    count = 1;
    firstSpikeBiPole = [];
    for opto = 1:size(optoTraceHit,1)
        optoId = find(diff(optoTraceHit(opto,:))==1)+1;
        for neuron = find(Spikes.BiPOLES.tagged==1)
            disp(['Calcuating neuron ', num2str(neuron) '...'])
            BiPolesNeuron = Spikes.PSTH.hit.spks{neuron};
            for optoIdx = 1:length(optoId)
                if optoIdx==length(optoId)
                    firstSpk = find(BiPolesNeuron(optoTrials(opto),optoId(optoIdx):end)==1,1);
                else
                    firstSpk = find(BiPolesNeuron(optoTrials(opto),optoId(optoIdx):(optoId(optoIdx+1)-10))==1,1);
                end
                if ~isempty(firstSpk)
                    firstSpikeBiPole(count) = firstSpk;
                    count = count+1;
                end
            end
        end
    end
    firstSpikeNoTagTotal{fileNum} = firstSpikeBiPole;
end
%% Plot it all out
combinedDat = horzcat(firstSpikeBiPoleTotal{:});
figure,
histogram(combinedDat,'binwidth',1,'normalization','probability'),axis square,set(gca,'tickdir','out'),box off
xlim([-5 50])
%% Fraction tagged
m_rt  = cellfun(@(x) sum(x)/length(x), taggedSpkTotal)';


col_noOpto = [0.4 0.4 0.4];        % gray base
col_pts_no = [0.2 0.2 0.2];        % dark gray points
col_pts_op = [0 123 167]/255;      % cerulean points [web:115][web:116]

figure; hold on

% Bars (semi‑transparent gray)
b = bar(1, mean(m_rt), 'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
b.FaceAlpha = 0.5;

% Error bars
% errorbar(1:2, m_rt, sem_rt, 'k', 'LineStyle','none', 'LineWidth',1);

% Overlay data points with horizontal jitter
jitter = 0.08;

x1 = 1 + (rand(size(m_rt))-0.5)*2*jitter;

scatter(x1, m_rt, 25, col_pts_no, 'filled', 'MarkerFaceAlpha',0.8);
% scatter(x2, optoHitrt,   25, col_pts_op, 'filled', 'MarkerFaceAlpha',0.8);

set(gca,'XTick',1,'XTickLabel',{'Opto'});
ylabel('Tagged Neurons (%)');
set(gca,'Box','off','TickDir','out','FontSize',12);
%%
neuronId  = 2;
tempSpk   = Spikes.PSTH.hit.spks{neuronId};   % [nTrials x nTime]
optoTrials = vertcat(IntanBehaviour.cueHitTrace.opto);   % [nTrials x 1], 0/1

nonOptoIdx = optoTrials == 0;
optoIdx = find(optoTrials == 1);

spk_nonOpto = tempSpk(nonOptoIdx,:);   % non‑opto trials

tempSpk   = Spikes.PSTH.hit.spks{neuronId};      % [nTrials x nTime], 0/1
optoTrials = vertcat(IntanBehaviour.cueHitTrace.opto);   % 0/1

spk_mod = tempSpk;        % copy to modify
optoAd      = 1;          % 0 = off, 1 = add opto-driven spikes
if optoAd==1
    spk_opto = adOpto(tempSpk,optoIdx);
else
    spk_opto    = tempSpk(optoIdx,:);      % opto trials
end
t = (1:size(tempSpk,2));               % time axis (samples or ms)

figure;

% -------- 1) Raster: non-opto --------
subplot(2,2,1); hold on
[row,col] = find(spk_nonOpto);
scatter(t(col), row, 6, [0.5 0.5 0.5], 'filled');                 % [web:6]
ylabel('Trials (no opto)');
title(['Example neuron: ', num2str(neuronId)]);

set(gca,'YDir','reverse','Box','off','TickDir','out');

% -------- 2) Raster: opto --------
subplot(2,2,2); hold on
[row,col] = find(spk_opto);
scatter(t(col), row, 6, [0 123 167] / 255, 'filled');         % magenta for opto
ylabel('Trials (opto)');
title('Optogenetic trials');
set(gca,'YDir','reverse','Box','off','TickDir','out');

% -------- 3) Mean rate: non-opto --------
binSize = 20;                           % samples or ms per bin
edges = 1:binSize:(size(tempSpk,2)+1);
centers = edges(1:end-1) + binSize/2;

% non-opto PSTH
cnt = zeros(size(spk_nonOpto,1), numel(edges)-1);
for b = 1:numel(edges)-1
    idx = edges(b):edges(b+1)-1;
    cnt(:,b) = sum(spk_nonOpto(:,idx),2);
end
rate_nonOpto = mean(cnt,1) * (1000/binSize);            % Hz [web:19]

subplot(2,2,3); hold on
plot(centers, rate_nonOpto, 'Color',[0.5 0.5 0.5],'LineWidth',1.5);
xlabel('Time (ms)');
ylabel('Rate (spks/s)');
title('Mean rate (no opto)');
set(gca,'Box','off','TickDir','out');

% -------- 4) Mean rate: opto --------
cnt = zeros(size(spk_opto,1), numel(edges)-1);
for b = 1:numel(edges)-1
    idx = edges(b):edges(b+1)-1;
    cnt(:,b) = sum(spk_opto(:,idx),2);
end
rate_opto = mean(cnt,1) * (1000/binSize);               % Hz [web:19]

subplot(2,2,4); hold on
plot(centers, rate_opto, 'Color',[0 123 167] / 255,'LineWidth',1.5);
xlabel('Time (ms)');
ylabel('Rate (spks/s)');
title('Mean rate (opto)');
set(gca,'Box','off','TickDir','out');
% Add vertical lines at each pulse
optoStart = 1500;      % ms
pulseStep = 50;        % ms between pulses
tEnd      = 3000;      % end of plotting window
pulseTimes = optoStart:pulseStep:tEnd;
for pt = pulseTimes
    xline(pt, '-', 'Color', [0.6 0.6 0.6]);   % dotted gray lines [web:68][web:88]
end

xlabel('Time (ms)');
ylabel('Rate (spks/s)');
title('Mean rate (opto)');
set(gca,'Box','off','TickDir','out');

%% LOCAL FUNCTIONS
function spk_opto = adOpto(tempSpk,optoIdx)
spk_mod = tempSpk;        % copy to modify


optoFreq    = 20;         % Hz
binSize_ms  = 1;          % your current spike bin size
optoStart   = 1500;       % first bin to consider (ms)
nPulses     = 20;          % how many pulses to simulate
stepBins    = round((1000/optoFreq)/binSize_ms);   % 50 bins

% Example: probability that can vary across pulses (bins)
% length must be >= nPulses
p_vec = linspace(0.9, 0.2, nPulses);   % low → high probability


for k = 1:numel(optoIdx)
    tr = optoIdx(k);

    pulse = 0;
    for t = optoStart:stepBins:size(spk_mod,2)
        pulse = pulse + 1;
        if pulse > numel(p_vec)
            break
        end

        p_this = p_vec(pulse);          % probability for this bin

        % draw spike for this bin in this trial
        if rand < p_this
            spk_mod(tr,t) = 1;
        end
    end
end


spk_opto    = spk_mod(optoIdx,:);      % opto trials
end

% helper to get post‑stim mean rate
function r = getRate(spk,binSize,stim_bin)
edges = 1:binSize:(size(spk,2)+1);
nBins = numel(edges)-1;
cnt = zeros(size(spk,1),nBins);
for b = 1:nBins
    idx = edges(b):edges(b+1)-1;
    cnt(:,b) = sum(spk(:,idx),2);
end
rate = mean(cnt,1) * (1000/binSize);
r = mean(rate(stim_bin:end));
end

