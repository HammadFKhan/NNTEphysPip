%% Batch Sequentiality Index Analysis

%% Parameters
parameters.Fs = 1000;
parameters.ts = 1/parameters.Fs;
parameters.windowBeforePull = 1.5;
parameters.windowAfterPull = 1.5;
parameters.windowBeforeCue = 1.5;
parameters.windowAfterCue = 1.5;
parameters.windowBeforeMI = 1.5;
parameters.windowAfterMI = 1.5;
parameters.experiment = 'cue';
parameters.cool = 0;
parameters.opto = 1;
parameters.eOPN = 1;
parameters.BiPOLES = 0;
NumEntropyBins = 20;

%% Select and read Excel file
% [excelName,excelFolder] = uigetfile({'*.xlsx;*.xls','Excel Files (*.xlsx, *.xls)'},'Select Excel file containing MAT file paths');
% if isequal(excelName,0), error('No Excel file selected.'); end
% excelFile = fullfile(excelFolder,excelName);

% M1 data during M1 inactivation
%excelFile = "\\10.165.57.13\Sutter_backup\Hammad\Ephys\LeverTask\Data_for_Figures\eOPN\M1Inactivation\eOPN_M1_Inactivation_Spikes.xlsx";

% M1 data during thalamic inactivation
excelFile = "\\10.165.57.13\Sutter_backup\Hammad\Ephys\LeverTask\Data_for_Figures\eOPN\ThalamusInactivation\eOPN_Th_Inactivation_Spikes.xlsx";

fileTable = readtable(excelFile,'VariableNamingRule','preserve','TextType','string');
nFiles = height(fileTable);
fprintf('\nFound %d files in Excel sheet.\n',nFiles);

%% Preallocate
SpikesBatch = repmat(struct('FileName',[],'FilePath',[],'animalID',[],'M1Spikes',[],'Hit',[],'nNeurons',[],'nTrials',[],'nTime',[],'Success',false,'Error',[]),nFiles,1);
FirstOptoTrialIndex = repmat(struct('FileName',[],'animalID',[],'Hit',NaN,'Miss',NaN,'FA',NaN),nFiles,1);
Seq = repmat(struct('FileName',[],'animalID',[],'baselineHit',[],'optoHit',[],'baselinePE',[],'optoPE',[],'baselineTS',[],'optoTS',[],'nBaselineHitTrials',[],'nOptoHitTrials',[],'Success',false,'Error',[]),nFiles,1);

%% Batch loop
for fileIdx = 1:nFiles

    fprintf('\n============================================================\n');
    fprintf('Processing file %d of %d\n',fileIdx,nFiles);
    fprintf('============================================================\n');

    fileName = fileTable.FileName(fileIdx);
    matFile = fileTable.FilePath(fileIdx);
    animalID = fileTable.animalID(fileIdx);

    SpikesBatch(fileIdx).FileName = char(fileName);
    SpikesBatch(fileIdx).FilePath = char(matFile);
    SpikesBatch(fileIdx).animalID = animalID;
    FirstOptoTrialIndex(fileIdx).FileName = char(fileName);
    FirstOptoTrialIndex(fileIdx).animalID = animalID;
    Seq(fileIdx).FileName = char(fileName);
    Seq(fileIdx).animalID = animalID;

    fprintf('File: %s | Animal: %s\n',fileName,animalID);

    %% Skip files without M1 spikes
    if fileTable.M1Spikes(fileIdx) == 0
        fprintf('Skipping: M1Spikes = 0\n');
        SpikesBatch(fileIdx).Error = 'M1Spikes = 0';
        Seq(fileIdx).Error = 'M1Spikes = 0';
        continue;
    end

    if ~isfile(matFile)
        warning('File not found: %s',matFile);
        SpikesBatch(fileIdx).Error = 'File not found';
        Seq(fileIdx).Error = 'File not found';
        continue;
    end

    try

        %% Check variables
        matVariables = who('-file',matFile);

        %% Load spikes
        if ismember('M1Spikes',matVariables)
            Dspikes = load(matFile,'M1Spikes');
            M1Spikes = Dspikes.M1Spikes;
            Spikes = M1Spikes;
        elseif ismember('Spikes',matVariables)
            Dspikes = load(matFile,'Spikes');
            Spikes = Dspikes.Spikes;
            M1Spikes = Spikes;
        else
            error('No M1Spikes or Spikes variable found.');
        end

%         SpikesBatch(fileIdx).M1Spikes = M1Spikes;

        %% Load/separate behaviour
        if ismember('IntanBehaviourBaseline',matVariables) && ismember('IntanBehaviourOpto',matVariables)
            Dbehaviour = load(matFile,'IntanBehaviourBaseline','IntanBehaviourOpto');
            IntanBehaviourBaseline = Dbehaviour.IntanBehaviourBaseline;
            IntanBehaviourOpto = Dbehaviour.IntanBehaviourOpto;
        elseif ismember('IntanBehaviour',matVariables)
            Dbehaviour = load(matFile,'IntanBehaviour');
            IntanBehaviour = Dbehaviour.IntanBehaviour;
            [IntanBehaviourBaseline,IntanBehaviourOpto,~,~] = separateOptoTrials(IntanBehaviour,parameters);
        else
            error('No IntanBehaviour variables found.');
        end

        %% First opto trial indices
        if isfield(IntanBehaviourBaseline,'cueHitTrace'), FirstOptoTrialIndex(fileIdx).Hit = numel(IntanBehaviourBaseline.cueHitTrace)+1; end
        if isfield(IntanBehaviourBaseline,'cueMissTrace'), FirstOptoTrialIndex(fileIdx).Miss = numel(IntanBehaviourBaseline.cueMissTrace)+1; end
        if isfield(IntanBehaviourBaseline,'MIFATrace'), FirstOptoTrialIndex(fileIdx).FA = numel(IntanBehaviourBaseline.MIFATrace)+1; end

        %% Check Hit spikes
        if ~isfield(Spikes,'PSTH') || ~isfield(Spikes.PSTH,'hit') || ~isfield(Spikes.PSTH.hit,'spks'), error('Spikes.PSTH.hit.spks not found.'); end
        if isempty(Spikes.PSTH.hit.spks), error('Spikes.PSTH.hit.spks is empty.'); end

        %% Convert spikes to trials x neurons x time
        nNeurons = numel(Spikes.PSTH.hit.spks);
        [nTrials,nTime] = size(Spikes.PSTH.hit.spks{1});
        hitSpikes = zeros(nTrials,nNeurons,nTime,'like',Spikes.PSTH.hit.spks{1});

        for nn = 1:nNeurons
            if ~isequal(size(Spikes.PSTH.hit.spks{nn}),[nTrials nTime]), error('Neuron %d has inconsistent Hit PSTH dimensions.',nn); end
            hitSpikes(:,nn,:) = Spikes.PSTH.hit.spks{nn};
        end

%         SpikesBatch(fileIdx).Hit = hitSpikes;
%         SpikesBatch(fileIdx).nNeurons = nNeurons;
%         SpikesBatch(fileIdx).nTrials = nTrials;
%         SpikesBatch(fileIdx).nTime = nTime;

        %% Split baseline and opto
        firstOptoHit = FirstOptoTrialIndex(fileIdx).Hit;
        if isnan(firstOptoHit), error('Could not determine first opto Hit trial.'); end
        if firstOptoHit < 2, error('No baseline Hit trials available.'); end
        if firstOptoHit > nTrials, error('No opto Hit trials available or trial counts do not match.'); end

        baselineHitSpikes = hitSpikes(1:firstOptoHit-1,:,:);
        optoHitSpikes = hitSpikes(firstOptoHit:end,:,:);

        baselineId = 1:length(IntanBehaviourBaseline.cueHitTrace);
        eOPNId = length(IntanBehaviourBaseline.cueHitTrace)+1:length(IntanBehaviour.cueHitTrace);
        assert(length(eOPNId)==length(IntanBehaviourOpto.cueHitTrace))
        % Since the sqEntropy is accessing each spike structure we need to seperate the trials
        baselineSpikes = Spikes;
        eOPNSpikes = Spikes;
        baselineSpikes.PSTH.hit.spks = cellfun(@(x) x(baselineId,:),baselineSpikes.PSTH.hit.spks,'UniformOutput',false);
        eOPNSpikes.PSTH.hit.spks = cellfun(@(x) x(eOPNId,:),eOPNSpikes.PSTH.hit.spks,'UniformOutput',false);

        %% Sequentiality Index
        % Rerun sequentiality based on optimal bins
        try
            Seq(fileIdx).baselineSqEntropy = getSqEntropy(baselineSpikes);
            Seq(fileIdx).eOPNSqEntropy = getSqEntropy(eOPNSpikes);
        catch ME
            disp('Error in looped Sequentiality ')

            [Seq(fileIdx).baselineHit,Seq(fileIdx).baselinePE,Seq(fileIdx).baselineTS] = SeqIndexDB(baselineHitSpikes,NumEntropyBins);
            [Seq(fileIdx).optoHit,Seq(fileIdx).optoPE,Seq(fileIdx).optoTS] = SeqIndexDB(optoHitSpikes,NumEntropyBins);
        end

        Seq(fileIdx).nBaselineHitTrials = size(baselineHitSpikes,1);
        Seq(fileIdx).nOptoHitTrials = size(optoHitSpikes,1);
        Seq(fileIdx).Success = true;
        SpikesBatch(fileIdx).Success = true;

%         fprintf('Baseline SI = %.4f | Opto SI = %.4f\n',Seq(fileIdx).baselineHit,Seq(fileIdx).optoHit);

    catch ME
        warning('Analysis failed for %s: %s',fileName,ME.message);
        SpikesBatch(fileIdx).Error = ME.message;
        Seq(fileIdx).Error = ME.message;
    end
end

% Summary
fprintf('\nBatch complete: %d/%d files successful.\n',sum([Seq.Success]),nFiles);


%% Plot Sequentiality Index: Baseline vs eOPN + Mixed Linear Effects Model

baseline = [];
opto = [];
animal = strings(0,1);
session = [];

for i = 1:numel(Seq)
    if Seq(i).Success && ~isempty(Seq(i).baselineHit) && ~isempty(Seq(i).optoHit)
        baseline(end+1,1) = Seq(i).baselineHit;
        opto(end+1,1) = Seq(i).optoHit;
        animal(end+1,1) = string(Seq(i).animalID);
        session(end+1,1) = i;
    end
end

% Mixed linear effects model
Y = [baseline;opto];
Condition = categorical([repmat("Baseline",numel(baseline),1);repmat("eOPN",numel(opto),1)],["Baseline","eOPN"]);
Animal = categorical([animal;animal]);
Session = categorical([session;session]);

tbl = table(Y,Condition,Animal,Session);

lme = fitlme(tbl,'Y ~ Condition + (1|Animal)');

disp(lme);
stats = anova(lme);
disp(stats);

coefTable = lme.Coefficients;
p = coefTable.pValue(strcmp(coefTable.Name,'Condition_eOPN'));

% Plot
figure('Color','w','Position',[100 100 240 320]); hold on;

x1 = 1;
x2 = 2;

for i = 1:numel(baseline)
    plot([x1 x2],[baseline(i) opto(i)],'-','Color',[0.65 0.65 0.65],'LineWidth',0.75);
end

scatter(x1*ones(size(baseline)),baseline,22,[0.7 0.7 0.7],'filled','MarkerEdgeColor','k','LineWidth',0.5);
scatter(x2*ones(size(opto)),opto,22,[1 0.35 0],'filled','MarkerEdgeColor','k','LineWidth',0.5);

plot(x1,mean(baseline),'k_','MarkerSize',16,'LineWidth',2);
plot(x2,mean(opto),'k_','MarkerSize',16,'LineWidth',2);

yl = ylim;
yrange = diff(yl);
yStat = yl(2)-0.04*yrange;

plot([x1 x1 x2 x2],[yStat-0.01*yrange yStat yStat yStat-0.01*yrange],'k','LineWidth',1);

if p < 0.0001
    pText = 'p < 0.0001';
else
    pText = sprintf('p = %.4f',p);
end

text(mean([x1 x2]),yStat+0.015*yrange,pText,'HorizontalAlignment','center','FontSize',9);

xlim([0.6 2.4]);
xticks([1 2]);
xticklabels({'Baseline','eOPN'});
ylabel('M1 Seq. Index');
box off;
set(gca,'TickDir','out','FontSize',9,'LineWidth',1);

draw now
%%
M1eOPN = Seq;
dat = [];
for n = 1:length(M1eOPN)
    try
        dat = vertcat(dat,[M1eOPN(n).baselineSqEntropy.CueHit.SqI(1)',M1eOPN(n).eOPNSqEntropy.CueHit.SqI(3)']);
    catch
        disp("No Sq detected")
        continue
    end
end
% dat = vertcat(dat,[M1eOPN(1).baselineSqEntropy.CueHit.SqI(3)',M1eOPN(4).eOPNSqEntropy.CueHit.SqI(6)']);
% [~, Id] = sort(diff(dat,1,2), 'ascend');
% dat = dat(Id(1:5),:);
plotDat(dat)
ylim([0.55 0.9])
axis square
%% Save
[saveName,saveFolder] = uiputfile('*.mat','Save Sequentiality Batch Results','M1SequentialityIndex_BatchResults_ThInactivation.mat');
if ~isequal(saveName,0)
    outputFile = fullfile(saveFolder,saveName);
    % save(outputFile,'SpikesBatch','FirstOptoTrialIndex','Seq','parameters','NumEntropyBins','fileTable','-v7.3');
    save(outputFile,'FirstOptoTrialIndex','Seq','parameters','NumEntropyBins','fileTable','-v7.3');
    fprintf('Saved: %s\n',outputFile);
end

function plotDat(dat)
baseline_data = dat(:,1);
eOPN_data = dat(:,2);
figure; % Create a new figure window

% Define colors for the points
baseline_color = [0.6 0.6 0.6]; % Gray
eOPN_color = [0.85 0.35 0.1];   % Orange

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
        'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', eOPN_color);

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
xticklabels({'Baseline', 'eOPN'}); % Set x-axis labels
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