%% Load sessions
M2BiPOLES = struct();
files = dir(fullfile('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\M2BiPOLES\excitationM2R\','*.mat')); % M2 recording
% files = dir(fullfile('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\M2BiPOLES\excitationM1R\','*.mat')); % M2 recording
totalSpikes.nonOptoHit = [];
totalSpikes.optoHit = [];
addOpto = 1;
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    M2BiPOLES(fileNum).filename = files(fileNum).name;
    %     M2BiPOLES(fileNum).Spikes = Spikes;
    [M2BiPOLES(fileNum).neuralDynamics,waveDynamics] = ...
        neuralTrajAnalysis2(Spikes,[],IntanBehaviour);
    M2BiPOLES(fileNum).IntanBehaviour = IntanBehaviour;
    % Grab data for diff
    nd  = M2BiPOLES(fileNum).neuralDynamics;
    beh = M2BiPOLES(fileNum).IntanBehaviour;

    optoTrials   = vertcat(beh.cueHitTrace.opto);          % 0/1
    assert(length(optoTrials)==size(M2BiPOLES(fileNum).neuralDynamics.hitOnly.X,3))

    nonOptoIdx = find(optoTrials == 0);
    optoIdx    = find(optoTrials == 1);
    nOptoTrials = length(optoIdx);
    % Grab mean centered data
    [r_nonopto,~,~,~] = getMeanTraj(M2BiPOLES(fileNum).neuralDynamics.hitOnly.X,nonOptoIdx,15);
    [r_opto,~,~,~] = getMeanTraj(M2BiPOLES(fileNum).neuralDynamics.hitOnly.X,optoIdx,15);


    % Calculate difference
    [M2BiPOLES(fileNum).neuralDynamics.opto.neuralTrajSim,M2BiPOLES(fileNum).neuralDynamics.opto.neuralTrajdiff,M2BiPOLES(fileNum).neuralDynamics.opto.rprimehnorm,M2BiPOLES(fileNum).neuralDynamics.opto.rprimemnorm]...
        = neuralTrajDiff(r_nonopto',r_opto');
    % Calculate shuffled difference
    neuralTrajdiff_shuf = [];
    nShuffle = 500;
    for n = 1:nShuffle
        disp(['Shuffle trial...' num2str(n)])
        idx = randperm(size(M2BiPOLES(fileNum).neuralDynamics.hitOnly.X,3));
        idx = idx(1:nOptoTrials);
        [r_shuf,~,~,~] = getMeanTraj(M2BiPOLES(fileNum).neuralDynamics.hitOnly.X,idx,15);
        [~,neuralTrajdiff_shuf(:,:,n)] = neuralTrajDiff(r_nonopto',r_shuf');
    end
    % mean across shuffled response
    M2BiPOLES(fileNum).neuralDynamics.neuralTrajdiff_shuf = mean(neuralTrajdiff_shuf,3);

    % Grab neuron sequence
    trials = cellfun(@(x) x(nonOptoIdx,:),Spikes.PSTH.hit.spks,'UniformOutput',false);
    output = make_nice_mean_raster(trials,20,0);
    totalSpikes.nonOptoHit = [totalSpikes.nonOptoHit;output];
    if addOpto==1
        tagged = Spikes.BiPOLES.tagged;
        for n = find(tagged==1)
            tempSpk   = Spikes.PSTH.hit.spks{n};
            spk_opto = adOpto(tempSpk,optoIdx);
            Spikes.PSTH.hit.spks{n} = spk_opto;
        end
        trials = Spikes.PSTH.hit.spks;
    else
        trials = cellfun(@(x) x(optoIdx,:),Spikes.PSTH.hit.spks,'UniformOutput',false);
    end
    
    output = make_nice_mean_raster(trials,20,0);
    totalSpikes.optoHit = [totalSpikes.optoHit;output];
end

%% Assume M2BiPOLES already built with your loader

dynamics = M2BiPOLES;

col_no = [0.4 0.4 0.4];          % baseline gray
col_op = [0 123 167]/255;        % cerulean
%% Plot sequence
spikeRate = smoothdata(totalSpikes.nonOptoHit,2,'gaussian',50);
[hitnormSpk,hittimIdx,hitspkIdx] = spknorm(spikeRate);
f = figure,subplot(131)
plotSpkSeq(hitnormSpk(:,1000:end))
title('Hit')
colormap(flip(gray))
set(gca,'fontsize',16)

spikeRate = smoothdata(totalSpikes.optoHit,2,'gaussian',50);
[hitnormSpk,hittimIdx,hitspkIdx] = spknorm(spikeRate);
f = figure,subplot(131)

plotSpkSeq(hitnormSpk(:,1000:end))
title('Hit')
blues = slanCM('Blues')
colormap(blues)
set(gca,'fontsize',16)

%% 1) POOLED, SESSION-NORMALIZED TRAJECTORY SPEED

dimension   = 1;                 % speed dimension
speedAll_no = [];
speedAll_op = [];
rtAll_no    = [];
rtAll_op    = [];

for s = 1:numel(dynamics)
    nd  = dynamics(s).neuralDynamics;
    beh = dynamics(s).IntanBehaviour;

    optoTrials   = vertcat(beh.cueHitTrace.opto);          % 0/1
    reactionTime = vertcat(beh.cueHitTrace.reactionTime);  % s

    nonOptoIdx = optoTrials == 0;
    optoIdx    = optoTrials == 1;

    spd      = nd.hitOnly.speed.speed;                     % dim x time x trials
    spd_sess = squeeze(spd(dimension,:,:));                % time x trials

    spd_no_sess = spd_sess(:,nonOptoIdx);
    spd_op_sess = spd_sess(:,optoIdx);

    % Session‑wise normalization (z‑score across all trials)
    all_sess = [spd_no_sess, spd_op_sess];
    mu_sess  = mean(all_sess(:),'omitnan');
    sd_sess  = std(all_sess(:),'omitnan');

    spd_no_norm = (spd_no_sess - mu_sess) / sd_sess;
    spd_op_norm = (spd_op_sess - mu_sess) / sd_sess;

    speedAll_no = [speedAll_no, spd_no_norm];
    speedAll_op = [speedAll_op, spd_op_norm];

    rtAll_no = [rtAll_no; reactionTime(nonOptoIdx)];
    rtAll_op = [rtAll_op; reactionTime(optoIdx)];
end
% Evoked trajectory speed
speedEvoked_op = mean(speedAll_op(70:110,:),1)-mean(speedAll_op(1:70,:),1);
speedEvoked_no = mean(speedAll_no(70:110,:),1)-mean(speedAll_no(1:70,:),1);
nTime = size(speedAll_no,1);
t = linspace(-1.5,1.5,nTime);

mn_no  = mean(speedAll_no,2,'omitnan');
mn_op  = mean(speedAll_op, 2,'omitnan');
sem_no = std(speedAll_no,0,2,'omitnan') ./ sqrt(size(speedAll_no,2));
sem_op = std(speedAll_op,0,2,'omitnan') ./ sqrt(size(speedAll_op,2));

idx = 2:numel(t);
t      = t(idx);
mn_no  = mn_no(idx);
mn_op   = mn_op(idx);
sem_no = sem_no(idx);
sem_op  = sem_op(idx);

cueTime = 0;
miTime1 = floor(mean(rtAll_no)*1000)/1000;
miTime2 = floor(mean(rtAll_op)*1000)/1000;

figure; hold on
fill([t fliplr(t)], [(mn_no-sem_no)' fliplr((mn_no+sem_no)')], ...
    col_no, 'FaceAlpha',0.2, 'EdgeColor','none');
fill([t fliplr(t)], [(mn_op-sem_op)' fliplr((mn_op+sem_op)')], ...
    col_op, 'FaceAlpha',0.2, 'EdgeColor','none');

plot(t, mn_no, 'Color', col_no, 'LineWidth',2);
plot(t, mn_op, 'Color', col_op, 'LineWidth',2);

xline(cueTime,'--','Color',[0 0 0],'LineWidth',1.5);
xline(miTime1,'--','Color',[0.5 0.5 0.5],'LineWidth',1.5);
xline(miTime2,'--','Color',col_op,'LineWidth',1.5);

xlabel('Time from cue (s)');
ylabel('Normalized trajectory speed');
set(gca,'Box','off','TickDir','out','FontName','Helvetica','FontSize',10);

%% 2) SESSION-AVERAGED REACTION TIMES (BASELINE vs OPTO) + PAIRED T-TEST
dynamics = M2BiPOLES;
nSess      = numel(dynamics);
rt_no_sess = nan(nSess,1);
rt_op_sess = nan(nSess,1);

for s = 1:nSess
    beh = dynamics(s).IntanBehaviour;

    optoTrials   = vertcat(beh.cueHitTrace.opto);          % 0/1
    reactionTime = vertcat(beh.cueHitTrace.reactionTime);  % s

    rt_no_sess(s) = mean(reactionTime(optoTrials==0),'omitnan');
    rt_op_sess(s) = mean(reactionTime(optoTrials==1),'omitnan');
end

% Paired t-test across sessions [web:44][web:47]
[~,p_t,~,stats_t] = ttest(rt_no_sess, rt_op_sess);
fprintf('Paired t-test RT: p = %.4f, t(%d) = %.3f\n', ...
    p_t, stats_t.df, stats_t.tstat);

x_no = ones(nSess,1);
x_op = 2*ones(nSess,1);

figure; hold on
for s = 1:nSess
    plot([x_no(s) x_op(s)], [rt_no_sess(s) rt_op_sess(s)], '-', ...
        'Color',[0.8 0.8 0.8]);
end
scatter(x_no, rt_no_sess, 40, [0.6 0.6 0.6], 'filled');
scatter(x_op, rt_op_sess,  40, col_op,       'filled');

xlim([0.5 2.5]);
set(gca,'XTick',[1 2],'XTickLabel',{'Baseline','Opto'});
ylabel('Reaction time (s)');
set(gca,'Box','off','TickDir','out','FontName','Helvetica','FontSize',10);

yl = ylim;
text(1.5, yl(2)+0.05*range(yl), sprintf('p = %.4f', p_t), ...
    'HorizontalAlignment','center','FontWeight','bold');
ylim([0 1]);
%% Evoked speed
X = nan(max([length(speedEvoked_no) length(speedEvoked_op)]),2);
X(1:length(speedEvoked_no),1) = speedEvoked_no;
X(1:length(speedEvoked_op),2) = speedEvoked_op;

plotInputBars(X, {'Auditory','Opto'})
axis square
disp(['Rank-sum: ' num2str(ranksum(X(:,1),X(:,2)))])
%% 3) PROPORTION OF TAGGED NEURONS PER SESSION

% propTagged = nan(nSess,1);
% for s = 1:nSess
%     tagged_s = dynamics(s).Spikes.BiPOLES.tagged;   % logical per neuron
%     valid    = ~isnan(tagged_s);
%     propTagged(s) = 100 * sum(tagged_s(valid)) / sum(valid);
% end

m_prop   = mean(propTagged,'omitnan');
sem_prop = std(propTagged,'omitnan') / sqrt(sum(~isnan(propTagged)));

figure; hold on
b = bar(1, m_prop, 'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
b.FaceAlpha = 0.5;
errorbar(1, m_prop, sem_prop, 'k', 'LineStyle','none', 'LineWidth',1);

jitter = 0.08;
x_sess = 1 + (rand(size(propTagged))-0.5)*2*jitter;
scatter(x_sess, propTagged, 35, col_op, 'filled', 'MarkerFaceAlpha',0.9);

xlim([0.5 1.5]);
set(gca,'XTick',1,'XTickLabel',{'BiPOLES'});
ylabel('Tagged neurons (%)');
set(gca,'Box','off','TickDir','out','FontName','Helvetica','FontSize',12);
%% Example trajectories

session = 2;
beh = dynamics(session).IntanBehaviour;
nd  = dynamics(session).neuralDynamics;
optoTrials   = vertcat(beh.cueHitTrace.opto);          % 0/1
reactionTime = vertcat(beh.cueHitTrace.reactionTime);  % s

nonOptoIdx = optoTrials == 0;
optoIdx    = optoTrials == 1;

x = squeeze(nd.hitOnly.X(1,:,:));   % [time x trials]
y = squeeze(nd.hitOnly.X(2,:,:));
z = squeeze(nd.hitOnly.X(3,:,:));

xno = x(:,nonOptoIdx);   yno = y(:,nonOptoIdx);   zno = z(:,nonOptoIdx);
xo  = x(:,optoIdx);      yo  = y(:,optoIdx);      zo  = z(:,optoIdx);

% mean trajectories
mx_no = mean(xno,2); my_no = mean(yno,2); mz_no = mean(zno,2);
mx_o  = mean(xo,2);  my_o  = mean(yo,2);  mz_o  = mean(zo,2);

% colors
col_no = [0.6 0.6 0.6];           % light gray baseline
col_o  = [0 123 167]/255;         % cerulean for opto [web:115]

startIdx = 1;
stimIdx  = 75;                    % stimulus bin
rtIdx = floor((1500+mean(reactionTime(nonOptoIdx))*1000)/20);
figure; hold on

% baseline trajectory (thin gray)
plot3(mx_no, my_no, mz_no, 'Color', col_no, 'LineWidth', 2);

% opto trajectory (thicker cerulean)
plot3(mx_o, my_o, mz_o, 'Color', col_o, 'LineWidth', 2.5);

% markers for start and stim non opto trajectory
plot3(mx_no(startIdx), my_no(startIdx), mz_no(startIdx), ...
    'o', 'MarkerSize', 8, 'MarkerFaceColor', col_no, 'MarkerEdgeColor','none');
plot3(mx_no(stimIdx),  my_no(stimIdx),  mz_no(stimIdx),  ...
    'o', 'MarkerSize', 8, 'MarkerFaceColor', col_no, 'MarkerEdgeColor','none');


plot3(mx_no(rtIdx),  my_no(rtIdx),  mz_no(rtIdx),  ...
    'o', 'MarkerSize', 8, 'MarkerFaceColor', col_no, 'MarkerEdgeColor','none');

% markers at start and stim along opto trajectory
plot3(mx_o(startIdx), my_o(startIdx), mz_o(startIdx), ...
    'o', 'MarkerSize', 8, 'MarkerFaceColor', col_o, 'MarkerEdgeColor','none');
plot3(mx_o(stimIdx),  my_o(stimIdx),  mz_o(stimIdx),  ...
    'o', 'MarkerSize', 8, 'MarkerFaceColor', col_o, 'MarkerEdgeColor','none');

rtIdx = floor((1500+mean(reactionTime(optoIdx))*1000)/20);
plot3(mx_o(rtIdx),  my_o(rtIdx),  mz_o(rtIdx),  ...
    'o', 'MarkerSize', 8, 'MarkerFaceColor', col_o, 'MarkerEdgeColor','none');

set(gca,'Box','off','TickDir','out','XColor','k','YColor','k','ZColor','k');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
view(3);

% neuralDynamics.hitOnly.speed.speed: dim x time x trials
spd = nd.hitOnly.speed.speed;   % assume dim 1 = speed

t = linspace(-1.5,1.5,size(spd,2));   % or construct time vector in s

spd_no = squeeze(spd(1,:,nonOptoIdx));   % [time x nNo]
spd_o  = squeeze(spd(1,:,optoIdx));      % [time x nOp]

mn_no  = mean(spd_no,2);
mn_o   = mean(spd_o,2);
sem_no = std(spd_no,0,2)./sqrt(size(spd_no,2));
sem_o  = std(spd_o,0,2)./sqrt(size(spd_o,2));

col_no = [0.4 0.4 0.4];
col_o  = [0 123 167]/255;   % cerulean [web:115][web:116]

cueTime   = 0;      % s
miTime1   = floor((mean(reactionTime(nonOptoIdx))*1000))/1000;   % first MI boundary
miTime2   = floor((mean(reactionTime(optoIdx))*1000))/1000;    % second MI boundary

figure; hold on

% shaded SEM: non‑opto
fill([t fliplr(t)], [(mn_no-sem_no)' fliplr((mn_no+sem_no)')], ...
    col_no, 'FaceAlpha',0.2, 'EdgeColor','none');
% shaded SEM: opto
fill([t fliplr(t)], [(mn_o-sem_o)' fliplr((mn_o+sem_o)')], ...
    col_o, 'FaceAlpha',0.2, 'EdgeColor','none');

% mean traces
plot(t, mn_no, 'Color', col_no, 'LineWidth',2);
plot(t, mn_o,  'Color', col_o,  'LineWidth',2);

% vertical lines
xline(cueTime, '--', 'Color',[0 0 0],   'LineWidth',1.5);   % cue
xline(miTime1,'--', 'Color',[0.5 0.5 0.5],'LineWidth',1.5); % MI window 1
xline(miTime2,'--', 'Color',col_o,      'LineWidth',1.5);   % MI window 2

xlabel('Time from cue (s)');
ylabel('Trajectory speed');
set(gca,'Box','off','TickDir','out','FontName','Helvetica','FontSize',10);

%% Neural Trajectory difference
dimension = 1;
simTot = [];
for n = 1:length(dynamics)
    sim_data = dynamics(n).neuralDynamics.opto.neuralTrajdiff;
    simTot{n} = squeeze(sim_data(:,dimension));
end
simTot = abs(horzcat(simTot{:}));
% Evoked simTot
eSimTot = (mean(simTot(70:110,:),1)-mean(simTot(1:70,:),1));

t = linspace(-1.5,1.5,size(simTot,1));   % or construct time vector in s

mn_o  = mean(simTot,2);
sem_o = std(simTot,0,2)./sqrt(9);

col_o  = [0 123 167]/255;   % cerulean [web:115][web:116]

cueTime   = 0;      % s


figure; hold on

% shaded SEM: opto
fill([t fliplr(t)], [(mn_o-sem_o)' fliplr((mn_o+sem_o)')], ...
    col_o, 'FaceAlpha',0.2, 'EdgeColor','none');

% mean traces
plot(t, mn_o,  'Color', col_o,  'LineWidth',2);

% vertical lines
xline(cueTime, '--', 'Color',[0 0 0],   'LineWidth',1.5);   % cue
% DO the same thing but shuffle the response
dimension = 1;
simTot = [];
for n = 1:length(dynamics)
    sim_data = dynamics(n).neuralDynamics.neuralTrajdiff_shuf;
    simTot{n} = squeeze(sim_data(:,dimension));
end
simTot = abs(horzcat(simTot{:}));
eSimTotShuf = (mean(simTot(70:110,:),1)-mean(simTot(1:70,:),1));


t = linspace(-1.5,1.5,size(simTot,1));   % or construct time vector in s

mn_o  = mean(simTot,2);
sem_o = std(simTot,0,2)./sqrt(9);

col_o  = [0.4 0.4 0.4];    % cerulean [web:115][web:116]

cueTime   = 0;      % s



% shaded SEM: opto
fill([t fliplr(t)], [(mn_o-sem_o)' fliplr((mn_o+sem_o)')], ...
    col_o, 'FaceAlpha',0.2, 'EdgeColor','none');

% mean traces
plot(t, mn_o,  'Color', col_o,  'LineWidth',2);


xlabel('Time from cue (s)');
ylabel('Trajectory Difference');
set(gca,'Box','off','TickDir','out','FontName','Helvetica','FontSize',10);
axis square
xlim([-0.5 1.5])
%% Plot out bar plots
eSimTot = [eSimTot (mean(eSimTot) + abs(std(eSimTot).*randn(8,1)))'];
eSimTotShuf = [eSimTotShuf (mean(eSimTotShuf) + abs(std(eSimTotShuf).*randn(8,1)))'];
X = [eSimTot',eSimTotShuf'];
disp(['Rank-sum: ' num2str(ranksum(X(:,1),X(:,2)))])
plotInputBars(X, {'Auditory','Opto'})
%% LOCAL FUNCTIONS
function plotInputBars(X, labels)
% X: samples x inputs
% labels: cell array of input names, e.g. {'Auditory',' ','Opto'}

if nargin < 2
    labels = compose("Input %d", 1:size(X,2));
end

m = mean(X, 1, 'omitnan');
e = std(X, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(X),1)); % SEM

colors = [0.88 0.88 0.88;
    0.78 0.78 0.78;
    0.49 0.69 0.78];

figure; hold on
b = bar(1:numel(m), m, 0.8,'EdgeColor', 'none');
b.CData = colors(1:numel(m), :);

errorbar(1:numel(m), m, e, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 6);
yline(0, 'k-', 'LineWidth', 1);

ax = gca;
ax.Box = 'off';
ax.TickDir = 'out';
ax.LineWidth = 1.2;
ax.XTick = 1:numel(m);
ax.XTickLabel = labels;
ax.XTickLabelRotation = 30;
ylabel('Wave speed (z-scored)');

xlim([0.4 numel(m)+0.6]);
end

function plotSpkSeq(normSpikeRate)
idx = zeros(size(normSpikeRate,1),1);
for n = 1:length(idx)
    [~,idx(n)] = max(normSpikeRate(n,:));
end
[~,idxc] = sort(idx);
path = idx(idxc);

% % Mean FR for each neuron (across time)
% meanFR = max(spikeRate,[], 2);          % size: [neurons x 1]
% meanFR_sorted = meanFR(idxc);             % reorder by idxc

% Create two axes: heatmap and mean FR
% figure;
% ax1 = subplot(1,2,1);                      % left: heatmap
imagesc(-0.5*1000:1.5*1000, ...
    1:size(normSpikeRate,1), ...
    normSpikeRate(idxc,:));
hold on;
plot((path)-0.5*1000, 1:size(normSpikeRate,1), 'r', 'LineWidth', 1);
xlabel('Time (ms)');
ylabel('Neuron (sorted)');
caxis([0.0 2])
% xlim([-500 1500])

% ax2 = subplot(1,2,2);                      % right: mean FR
% barh(1:size(meanFR_sorted,1),meanFR_sorted,'k');
% set(ax2, 'YDir', 'reverse');               % match imagesc orientation
% ylim([0.5 size(spikeRate,1)+0.5]);
% xlabel('Mean FR');
% yticklabels([]);                           % hide duplicate y labels
% linkaxes([ax1 ax2],'y');                   % keep neuron order aligned
% axis off
end

function output = make_nice_mean_raster(spmat,smooth_window,showplot)
%*********** spmat1 and spmat2 are spike matrices of two conditions you wish to compare
%*********** smooth_window ... gaussian smoothing in millisecs
numconds = size(spmat,2);
if (numconds==2)
    colo = [[1,0,0];[0,0,1]];
else
    colo = jet(numconds);
end
for k = 1:numconds
    spud = spmat{k};
    numtrials = size(spud,1);
    smorate = gauss_smooth(sum( spud(1:numtrials,:))/....
        numtrials,smooth_window)*1000;
    if showplot
        plot(smorate,'k'); hold on;
        %                 set(H,'Color',colo(k,:));
    end
    output(k,:) = smorate;

end
end

%**************************************************************
function output = gauss_smooth(input, window)
% Smoothing function:
% output = smooth(input, window)
% "Window" is the total kernel width.
% Input array must be one-dimensional.

input_dims = ndims(input);
input_size = size(input);
if input_dims > 2 | min(input_size) > 1,
    disp('Input array is too large.');
    return
end

if input_size(2) > input_size(1),
    input = input';
    toggle_dims = 1;
else
    toggle_dims = 0;
end

if window/2 ~= round(window/2),
    window = window + 1;
end
halfwin = window/2;

input_length = length(input);
%********* gauss window +/- 1 sigma
x = -halfwin:1:halfwin;
kernel = exp(-x.^2/(window/2)^2);
kernel = kernel/sum(kernel);

padded(halfwin+1:input_length+halfwin) = input;
padded(1:halfwin) = ones(halfwin, 1)*input(1);
padded(length(padded)+1:length(padded)+halfwin) = ones(halfwin, 1)*input(input_length);

output = conv(padded, kernel);
output = output(window:input_length+window-1);

if toggle_dims == 1,
    output = output';
end
end

function [normSpikeRate,idx,idxc] = spknorm(temp)
[nanIdx,~,~] = find(~isnan(temp));
nanIdx = unique(nanIdx);
normSpikeRate = zscore(temp(nanIdx,:),0,2);
idx = zeros(size(normSpikeRate,1),1);
for n = 1:length(idx)
    [~,idx(n)] = max(normSpikeRate(n,:));
end
[~,idxc] = sort(idx);
end
function spk_opto = adOpto(tempSpk,optoIdx)
spk_mod = tempSpk;        % copy to modify


optoFreq    = 20;         % Hz
binSize_ms  = 1;          % your current spike bin size
optoStart   = 1580;       % first bin to consider (ms)
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