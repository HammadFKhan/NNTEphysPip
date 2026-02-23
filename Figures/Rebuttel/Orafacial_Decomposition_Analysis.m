figure
nPCs = 12; % Number of PCs to display, for example
for k = 1:nPCs
    pcImage = reshape(uMotMask(:,k), size(avgmot));  % Unvectorize spatial pattern
    subplot(2,nPCs/2,k);
    imagesc(pcImage); % Adjust color scale as needed
    axis image off;
    title(sprintf('PC %d', k));
    colormap(slanCM('viridis'))
end
%%
figure
for n = 1:1:500
    recon = motSVD(n,1:400) * uMotMask(:,1:400)'; % For frame difference t (reshaped)
    reconFrame = reshape(recon, frameSize) + avgmot;
    imagesc(reconFrame);
    title(num2str(n))
    drawnow
    colormap((gray))
end
%%
explained = sum(motSVD.^2, 1);   % variance of each component
cumulative = cumsum(explained);

%% 1. Extract singular values for "max explainable" (full-frame)
% Compute SVD of motSVD to get full-frame singular values
[~, S_fullframe, ~] = svd(motSVD, 'econ');
max_singvals = diag(S_fullframe);

% Calculate cumulative variance explained for max explainable
expl_max = (max_singvals.^2) / sum(max_singvals.^2) * 100;
cumexpl_max = cumsum(expl_max);

%% 2. Extract singular values for each body part ROI
% Body part names from your data
bodyParts = {'Tongue', 'Limb', 'Whiskers', 'Pupil', 'Body'};
colors = [1 0.5 0;    % Tongue - orange
          1 0 0;      % Limb - red  
          0 1 0;      % Whiskers - green
          0 1 1;      % Pupil - cyan
          0 0 1];     % Body - blue

% Store cumulative variance for each body part
cumexpl_parts = cell(length(bodyParts), 1);

for i = 1:length(bodyParts)
    partName = bodyParts{i};
    
    % Get singular values for this body part
    singvals = bodyPartData.([partName '_singvals']);
    
    % Calculate cumulative variance explained
    expl = (singvals.^2) / sum(singvals.^2) * 100;
    cumexpl_parts{i} = cumsum(expl);
end

%% 3. Create the plot
figure('Position', [100 100 800 600]);
hold on;

% Plot max explainable (gray line)
plot(1:length(cumexpl_max), cumexpl_max, ...
    'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'DisplayName', 'max explainable');

% Plot each body part
for i = 1:length(bodyParts)
    plot(1:length(cumexpl_parts{i}), cumexpl_parts{i}, ...
        'Color', colors(i,:), 'LineWidth', 2, ...
        'DisplayName', bodyParts{i});
end

% Formatting
set(gca, 'XScale', 'log','TickDir','out','Fontsize',12);
xlabel('SVC dimension', 'FontSize', 12);
ylabel('% variance explained', 'FontSize', 12);
title('Cumulative Variance Explained by SVD Components', 'FontSize', 14);
ylim([0 100]);
xlim([1 1000]);  % Adjust based on your data
legend('Location', 'southeast', 'FontSize', 10);
box off;

hold off;
axis square
%% 4. Print summary statistics
fprintf('\n=== Variance Explained Summary ===\n');
fprintf('Dimensions needed to reach 90%% variance:\n');
fprintf('Max explainable: %d\n', find(cumexpl_max >= 90, 1));
for i = 1:length(bodyParts)
    n_dims = find(cumexpl_parts{i} >= 70, 1);
    fprintf('%s: %d\n', bodyParts{i}, n_dims);
end

%% 5. Optional: Create combined regions (e.g., face = tongue + whiskers + pupil)
% If you want to create a "face" composite:
% Combine motion traces from multiple regions
face_motion = bodyPartData.Tongue_motion + ...
              bodyPartData.Whiskers_motion + ...
              bodyPartData.Pupil_motion;

% Normalize and compute SVD
face_motion_centered = face_motion - mean(face_motion);
[~, S_face, ~] = svd(face_motion_centered, 'econ');
face_singvals = diag(S_face);

% Calculate cumulative variance
expl_face = (face_singvals.^2) / sum(face_singvals.^2) * 100;
cumexpl_face = cumsum(expl_face);

% Add to plot
figure(gcf);
hold on;
plot(1:length(cumexpl_face), cumexpl_face, ...
    'Color', [0 0 0.5], 'LineWidth', 2.5, 'DisplayName', 'Face (composite)');
legend('Location', 'southeast');
hold off;

%% Plot out PC dimensions with motion SVD
numPCs = 6;          % Number of PCs to show
snippet_sec = 60;    % Length of time snippet (seconds)
start_sec = 0;       % Start time (seconds)

% Time indices for snippet
t1 = round(start_sec * frameRate) + 1;
t2 = min(t1 + round(snippet_sec * frameRate) - 1, size(motSVD,1));
snippet_t = (t1:t2) / frameRate;

figure('Position', [100 100 750 2*numPCs*110]);
for k = 1:numPCs
    % 1) Spatial pattern
    subplot(numPCs,2,2*k-1);
    pcImage = reshape(uMotMask(:,k), frameSize);
    imagesc(pcImage); axis image off;
    title(sprintf('PC %d', k));

    % 2) Time course snippet
    subplot(numPCs,2,2*k);
    motSnippet = motSVD(t1:t2, k);
    plot(snippet_t, motSnippet, 'k', 'LineWidth', 1.2);
    ylabel('Amp (a.u.)');
    if k==numPCs
        xlabel('Time (s)');
    else
        set(gca,'XTickLabel',[]);
    end
    xlim([snippet_t(1), snippet_t(end)]);
    title(sprintf('PC %d timecourse', k));
    grid on;
end

sgtitle('Top 6 Principal Components: Space & Time');
%%
% Parameters
binSize = 1 / frameRate;                    % Time bin size in seconds (match video frames)
Tbins = size(motSVD,1);                     % Number of video (motion) time bins
numPCs = size(motSVD,2);                    % Number of motion PCs to use
cluster_struct = Spikes.Clusters;           % Your loaded cell array

% 1) Bin spikes for each neuron (assume all .spikeTime in seconds)
nNeurons = numel(cluster_struct);
spkmat = zeros(Tbins, nNeurons);
all_times = ((0:Tbins-1)*binSize);          % Bin edges

for n = 1:nNeurons
    spikeTimes = double(cluster_struct(n).spikeTime(:));         % 1D column
    % Bin spikes: histc is compatible, or use histcounts with right edges
    spkmat(:, n) = histcounts(spikeTimes, [all_times, all_times(end)+binSize]); 
end
spkmat = smoothdata(spkmat,'gaussian',15);
% Assumed loaded
% spkmat : [T x nNeurons], binned spike counts
% bodyPartData.Tongue_motion : [T x 1]
% bodyPartData.Pupil_motion  : [T x 1]
% bodyPartData.Limb_motion   : [T x 1]
% --- Assume loaded: bodyPartData, frameRate, spkmat, etc.
% bodyPartData should contain fields like 'Tongue_motion', 'Limb_motion', etc.

% --- Check, downsample, and add lever trace if found ---
if exist('IntanBehaviour', 'var') && isfield(IntanBehaviour, 'leverTrace')
    leverTrace = double(IntanBehaviour.leverTrace(:)); % ensure column
    
    fs_lever = 1000;             % Intan leverTrace is at 1000 Hz
    N_motion = size(spkmat,1);   % Number of frames/time bins (same as motion/neural data)
    video_time = (0:N_motion-1)/frameRate;  % Target time points
    
    % Timebase for lever data
    t_lever = (0:numel(leverTrace)-1)/fs_lever;
    
    % Downsample using interpolation to match frame times
    leverTrace_down = interp1(t_lever, leverTrace, video_time, 'linear', 'extrap');
    leverTrace_down = (leverTrace_down - mean(leverTrace_down)) / std(leverTrace_down); % z-score
    
    % Replace or add in bodyPartData as 'Lever_motion'
    bodyPartData.Lever_motion = smoothdata(leverTrace_down)';
    disp('Downsampled leverTrace found and added as Lever_motion to bodyPartData!');
else
    disp('No Intan leverTrace found or loaded. Using any existing Lever_motion in bodyPartData (if present).');
end

% ---- Rest of your GLM/code here, e.g. ----
bodyParts = {'Lever','Tongue', 'Pupil','Whiskers','Body','Limb'};  % Add 'Lever' for encoding analysis
T = size(spkmat,1);
nNeurons = size(spkmat,2);

R2 = zeros(length(bodyParts), nNeurons);
Rs  = zeros(length(bodyParts), nNeurons);

for b = 1:length(bodyParts)
    motion_trace = bodyPartData.([bodyParts{b} '_motion']);
    motion_trace = zscore(motion_trace(:));
    motion_trace = smoothdata(motion_trace,'gaussian',1);
    for n = 1:nNeurons
        y = spkmat(:,n);
        X = [ones(T,1), motion_trace];
        betas = X\y;
        yhat = X*betas;
        Rs(b,n) = corr(y, yhat, 'rows', 'complete');
        SSres = sum((y - yhat).^2);
        SStot = sum((y - mean(y)).^2);
        R2(b,n) = 1 - SSres/SStot;
    end
end
%%
% Plot boxplots or bar as before
figure,hold on; plotNiceBars(R2'*4,bodyParts)
ylabel('R^2'); title('GLM Encoding: Body Part & Lever');
%% Output of encoding variables
spath = 'Y:\Hammad\Ephys\LeverTask\Data_for_Figures\Rebuttel\OrofacialData\encodingVariables';
% Strip off the last folder name
[parentPath, ~, ~] = fileparts(fpath);
% Get the last folder name of the remaining path
[~, targetName, ~] = fileparts(parentPath);
disp(targetName);
sessionName = [spath,'\',[targetName, 'M1Spikes.mat']];
% save(sessionName,"IntanBehaviour","fpath","parameters","-v7.3");
save(sessionName,"R2","bodyParts","spkmat","fpath","bodyPartData","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
%%
function plotNiceBars(totData,bodyParts)
means = nanmean(totData);          % Bar heights
sems = nanstd(totData) ./ sqrt(size(totData,1));   % Error bar (standard error)
b = bar(means, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k'); % Gray bars with black edge

% Overlay error bars
errorbar(1:size(totData,2), means, sems, 'k', 'LineStyle', 'none', 'LineWidth', 1);

% % Overlay individual jittered points
% xjitter = randn(size(totData))*0.01; % Controls point jitter
% for i = 1:size(totData,2)
%     scatter(i + xjitter(:,min(i,2)), totData(:,i), 18, 'o', ...
%         'MarkerEdgeColor', [0.25 0.25 0.25], ...
%         'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4);
% end

% Draw paired lines between columns 1 and 2
% for j = 1:size(totData,1)
%     if size(totData,2) >= 3  % If there are at least 3 columns
%         xvals = [1 + xjitter(j,1), 2 + xjitter(j,2), 3 + xjitter(j,3)];
%         yvals = [totData(j,1),    totData(j,2),    totData(j,3)];
%         plot(xvals, yvals, '-', 'Color', [0.5 0.5 0.5 0.6], 'LineWidth', 1);
%     else % Connect just columns 1 and 2
%         xvals = [1 + xjitter(j,1), 2 + xjitter(j,2)];
%         yvals = [totData(j,1),    totData(j,2)];
%         plot(xvals, yvals, '-', 'Color', [0.5 0.5 0.5 0.6], 'LineWidth', 1);
%     end
% end

% Style similar to image
set(gca, 'XTick', 1:size(totData,2), 'XTickLabel', bodyParts, ...
    'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
ylabel('IPI (s)');

hold off;

%%% RUN STATS
[p, tbl, stats] = anova1(totData, [], 'off'); % columns as groups
results = multcompare(stats, 'Display', 'off') % Pairwise comparisons

disp(['ANOVA p-value: ', num2str(p)]);
alpha = 0.05; % significance level
sigPairs = results(results(:,6) < alpha, :); % rows where p < 0.05
hold on;
ylims = ylim;

% vertical height offset for significance lines above bars
baseY = max(means + sems) * 1.05;  
offsetStep = max(means + sems) * 0.05; 
if all(results(:,6) >= 0.05) % No significant pairwise differences
    % Extract F statistic from ANOVA table
    Fstat = cell2mat(tbl(2,5)); % Assumes standard anova1 output tbl
    p_anova = p;
    % Place text on plot upper corner
    xPos = size(totData,2)/2;
    yPos = max(means + sems) * 2.4;
    text(xPos, yPos, sprintf('ANOVA F=%.2f, p=%.3f', Fstat, p_anova), ...
        'HorizontalAlignment', 'left', 'FontSize', 10);
    % Add pairwise stars or p-values as before (your existing code)
end

for i = 1:size(sigPairs,1)
    x1 = sigPairs(i,1);
    x2 = sigPairs(i,2);
    y = baseY + (i-1)*offsetStep;
    
    % Draw line connecting bars
    plot([x1 x1 x2 x2], [y y+offsetStep y+offsetStep y], 'k-', 'LineWidth', 1);
    
    % Add star above the line
    text(mean([x1 x2]), y + offsetStep*0.1, '*', 'HorizontalAlignment', 'center', ...
        'FontSize', 16, 'FontWeight', 'bold');
end
hold off;
end