%% BiPOLES Sq index 
M1baseline = struct();
M1opto = struct();
% files = dir(fullfile('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\M2BiPOLES\excitationM1R\','*.mat')); % M1 Recording
% files = vertcat(files,dir(fullfile('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\M2BiPOLES\excitationM2R\','*.mat'))); % M2 recording
files = dir(fullfile('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\M2BiPOLES\excitationM2R\','*.mat')); % M1 Recording

for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    % Now calculate time from first spike based on opto tag M1 neurons and the
    % time in which the pulse reaches one.
    % BiPOLES index is based off the hit PSTH
    optoId = vertcat(IntanBehaviour.cueHitTrace.opto);   % 0/1 labels
    optoTrials = find(optoId==1);
    nonoptoTrials = find(optoId==0);
    % Do some assertions so we know we got it right
    assert(size(vertcat(IntanBehaviour.cueHitTrace.opto),1)==size(Spikes.PSTH.hit.spks{1},1))
    BiPolesNeuron = struct();
    count = 1;
    BiPolesNeuron.PSTH = Spikes.PSTH;
    BiPolesNeuron.PSTH.hit.spks = [];
    for neuron = find(Spikes.BiPOLES.tagged==1)
        disp(['Calcuating neuron ', num2str(neuron) '...'])
        BiPolesNeuron.PSTH.hit.spks{count} = Spikes.PSTH.hit.spks{neuron}(optoTrials,:);
        count = count+1;
    end
    % Grab only BiPole trials and non tagged neurons
    BiPoles = struct();
    BiPoles.PSTH = Spikes.PSTH;
    BiPoles.PSTH.hit.spks = [];
    for neuron = 1:length(Spikes.PSTH.hit.spks)
        disp(['Calcuating neuron ', num2str(neuron) '...'])
        BiPoles.PSTH.hit.spks{neuron} = Spikes.PSTH.hit.spks{neuron}(optoTrials,:);
    end
    % now grab the non tagged
    count = 1;
    nonBiPolesNeuron = struct();
    nonBiPolesNeuron.PSTH = Spikes.PSTH;
    nonBiPolesNeuron.PSTH.hit.spks = [];
    for neuron = find(Spikes.BiPOLES.tagged==0)
        disp(['Calcuating neuron ', num2str(neuron) '...'])
        nonBiPolesNeuron.PSTH.hit.spks{count} = Spikes.PSTH.hit.spks{neuron}(nonoptoTrials,:);
        count = count+1;
    end
    % Grab all neurons on nonopto trial
    nonBiPoles = struct();
    nonBiPoles.PSTH = Spikes.PSTH;
    nonBiPoles.PSTH.hit.spks = [];
    for neuron = 1:length(Spikes.PSTH.hit.spks)
        disp(['Calcuating neuron ', num2str(neuron) '...'])
        nonBiPoles.PSTH.hit.spks{neuron} = Spikes.PSTH.hit.spks{neuron}(nonoptoTrials,:);
    end



    M1baseline(fileNum).allNeuron.baselineSqEntropy = getSqEntropy(nonBiPoles);
    M1baseline(fileNum).taggedNeuron.baselineSqEntropy = getSqEntropy(nonBiPolesNeuron);

    M1opto(fileNum).allNeuron.optoSqEntropy = getSqEntropy(BiPoles);
    M1opto(fileNum).taggedNeuron.optoSqEntropy = getSqEntropy(BiPolesNeuron);
end
%% Plot it out
dat = [];
for n = 1:length(M1baseline)
    if isfield(M1baseline(n).allNeuron.baselineSqEntropy.CueHit,'SqI') && isfield(M1opto(n).taggedNeuron.optoSqEntropy.CueHit,'SqI')
        dat = vertcat(dat,[M1baseline(n).allNeuron.baselineSqEntropy.CueHit.SqI(2)',M1opto(n).taggedNeuron.optoSqEntropy.CueHit.SqI(2)']);
        dat = vertcat(dat,[M1baseline(n).allNeuron.baselineSqEntropy.CueHit.SqI(3)',M1opto(n).taggedNeuron.optoSqEntropy.CueHit.SqI(3)']);
    end
end
% [~, Id] = sort(diff(dat,1,2), 'ascend');
% dat = dat(Id(1:5),:);
plotDat(dat)
ylim([0.4 0.9])
axis square

%%
dat = [];
for n = 1:length(M1baseline)
    if isfield(M1baseline(n).baselineSqEntropy.CueHit,'PE') && isfield(M1cooling(n).cooledSqEntropy.CueHit,'PE')
    dat = vertcat(dat,[M1baseline(n).baselineSqEntropy.CueHit.PE(1)',M1cooling(n).cooledSqEntropy.CueHit.PE(2)']);
    end
end
% [~, Id] = sort(diff(dat,1,2), 'ascend');
% dat = dat(Id(1:5),:);
plotDat(dat)
ylim([0.4 0.9])
axis square
ylabel('PE', 'FontSize', 16); 


dat = [];
for n = 1:length(M1baseline)
    if isfield(M1baseline(n).baselineSqEntropy.CueHit,'TS') && isfield(M1cooling(n).cooledSqEntropy.CueHit,'TS')
    dat = vertcat(dat,[M1baseline(n).baselineSqEntropy.CueHit.TS(1)',M1cooling(n).cooledSqEntropy.CueHit.TS(2)']);
    end
end
% [~, Id] = sort(diff(dat,1,2), 'ascend');
% dat = dat(Id(1:5),:);
plotDat(dat)
ylim([0.8 1])
axis square
ylabel('TS', 'FontSize', 16); 

%%
function plotDat(dat)
baseline_data = dat(:,1);
eOPN_data = dat(:,2);
figure; % Create a new figure window

% Define colors for the points
baseline_color = [0.6 0.6 0.6]; % Gray
cooling_color = [0.1 0.35 0.85];   % Orange

% Plot lines connecting paired points first (light gray, thin)
hold on; % Keep the plot active for multiple elements
for i = 1:numel(baseline_data)
    plot([1, 2], [baseline_data(i), eOPN_data(i)], 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5); % Light gray line
end

% Plot individual data points using scatter
% Baseline points (x=1)
scatter(ones(size(baseline_data)), baseline_data, 100, 'filled', ...
        'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', baseline_color);

% eOPN points (x=2)
scatter(2 * ones(size(eOPN_data)), eOPN_data, 100, 'filled', ...
        'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', cooling_color);

% Add a horizontal dashed line at y = 0 for reference
%plot(xlim, [0 0], 'k--', 'LineWidth', 1.5);

% --- Beautify the plot ---
ax = gca; % Get current axes handle

ax.TickDir = 'out'; % Ticks point outwards
ax.FontSize = 14;   % Font size for tick labels
ax.Box = 'off';     % Turn off the box around the plot

% Set x-axis limits and labels for two groups
xlim([0.5 2.5]); % Adjust limits to center the two columns
xticks([1 2]); % Set tick marks at x=1 and x=2
xticklabels({'Nonopto', 'Opto'}); % Set x-axis labels
ax.XAxis.Color = [0.3 0.3 0.3]; % Darker color for x-axis labels
ax.YAxis.Color = [0.3 0.3 0.3]; % Darker color for y-axis labels
xlabel(''); % No overall x-axis label needed

% Y-axis label (match the image more closely)
ylabel('SI Index', 'FontSize', 16); 

ylim([0.6 0.9]); % Adjust limits to ensure 0 is visible

% --- Add statistical annotation (line and p-value) ---
% Get current y-axis limits to place the p-value
yLimits = ylim(ax);
xLimits = xlim(ax);

% Position for the p-value line and text (adjust these values manually for best fit)
y_line = yLimits(2) * 0.95; % Near the top
x_left = 1; % Corresponds to Baseline x-position
x_right = 2; % Corresponds to eOPN x-position

% Draw the line connecting the two groups for annotation
line([x_left, x_right], [y_line, y_line], 'Color', 'k', 'LineWidth', 0.5);

% Draw the small vertical bars at the ends of the horizontal line
line([x_left, x_left], [y_line - (range(yLimits)*0.02), y_line], 'Color', 'k', 'LineWidth', 0.5);
line([x_right, x_right], [y_line - (range(yLimits)*0.02), y_line], 'Color', 'k', 'LineWidth', 0.5);

% Add the p-value text (using the paired t-test p-value as in your second image)
text_x_pos = (x_left + x_right) / 2;
text_y_pos = y_line + (range(yLimits)*0.03); % Slightly above the line
[h_ttest_paired, p_ttest_paired, ci_ttest_paired, stats_ttest_paired] = ttest2(eOPN_data, baseline_data);
text(text_x_pos, text_y_pos, sprintf('p = %.4f', p_ttest_paired), ...
     'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'bottom', ...
     'FontSize', 14, 'FontWeight', 'bold'); % Similar font size/weight as image

% Optional: Add a title if you want, but the image you provided doesn't have one
% title('TW Speed Modulation', 'FontSize', 16);

hold off; % Release the plot
end
