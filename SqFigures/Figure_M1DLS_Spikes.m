% Warp response pooled together for M1 and DLS
M1Sq= struct();
DLSSq = struct();
files = dir(fullfile('D:\SequenceProject\WarpedSpikes\M1\','*.mat'));
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    M1Sq(fileNum).warpedSpikes = warpedSpks;
    M1Sq(fileNum).IntanBehaviour = IntanBehaviour;
    M1Sq(fileNum).filename = files(fileNum).name;
end

files = dir(fullfile('D:\SequenceProject\WarpedSpikes\DLS\','*.mat'));
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    DLSSq(fileNum).warpedSpikes = warpedSpks;
    DLSSq(fileNum).IntanBehaviour = IntanBehaviour;
    DLSSq(fileNum).filename = files(fileNum).name;
end
%% Get statistics
dynamics = DLSSq;
modulationIndex = [];
selectivity_index = [];
best_pull = [];
for n = 1:length(dynamics)
dynamics(n).warpedSpikes = getWarpSpkStats(dynamics(n).warpedSpikes);
modulationIndex = vertcat(modulationIndex,dynamics(n).warpedSpikes.stats.modulationIndex);
selectivity_index = vertcat(selectivity_index,dynamics(n).warpedSpikes.stats.selectivityIndex);
best_pull = vertcat(best_pull,dynamics(n).warpedSpikes.stats.bestPull);
end
num_pulls = size(modulationIndex,2);
%% Plot Modulation for each pull

bin_edges = -1:0.1:1; % Adjust bin size if desired
colors = [0 0 1; 0 0.5 0; 0.7 0 0]; % blue, green, red for pulls
figure; hold on;

for p = 1:num_pulls
    subplot(1,3,p),
    % Bin and normalize
    counts = histcounts(modulationIndex(:,p), bin_edges, 'Normalization', 'probability');
    bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
    % Stairs plot
    stairs(bin_centers, counts, 'Color', colors(p,:), 'LineWidth', 2);
    xlabel('Modulation Index');
    ylabel('Probability');
    title('Modulation Index');
    hold off;
    set(gca,'tickdir','out'),axis square
    ylim([0 0.2])
    box off
end


%% temp(1:length(sel_for_pull),p) = sel_for_pull;
temp = nan(length(selectivity_index),3);
y = [];
for p = 1:num_pulls
    sel_for_pull = selectivity_index(best_pull == p); % Example: your selectivity indices for this pull
    temp(1:length(sel_for_pull),p) = sel_for_pull;
    y(p) = length(sel_for_pull)/length(selectivity_index);
end

color = [46,49,179;46,149,49;179,49,46]/255;
figure,violinplot(temp,[],'ViolinColor',color);
xlabel('Pull Preference (Selectivity) Index')
ylabel('Number of Neurons')
legend({'Pull 1', 'Pull 2', 'Pull 3'}, 'Location', 'Best')
title('Distribution of Neuronal Pull Preference')
hold off;box off,set(gca,'tickdir','out'),axis square
figure,
x = [1];
bar(x,y,'stacked')
hold off;box off,set(gca,'tickdir','out'),axis square
xlim([0 2])
%%
% 
figure;
subplot(131),scatter(responseTot(:,1), responseTot(:,3), 24, 'filled','k');
hold on;
plot([0 1], [0 1], 'k--', 'LineWidth', 1.2); % y = x reference line
hold off;
xlabel('Mean Spike Response', 'FontWeight', 'bold');
ylabel('Mean Spike Response', 'FontWeight', 'bold');
title('Pull 1 vs Pull 3', 'FontWeight', 'bold');
set(gca, 'FontSize', 8, 'Box', 'off', 'GridAlpha', 0.4);
axis square;
box off,set(gca,'tickdir','out')

subplot(132),scatter(responseTot(:,1), responseTot(:,2), 24, 'filled','k');
hold on;
plot([0 1], [0 1], 'k--', 'LineWidth', 1.2); % y = x reference line
hold off;
xlabel('Mean Spike Response', 'FontWeight', 'bold');
ylabel('Mean Spike Response', 'FontWeight', 'bold');
title('Pull 1 vs Pull 2', 'FontWeight', 'bold');
set(gca, 'FontSize', 8, 'Box', 'off', 'GridAlpha', 0.4);
axis square;
box off,set(gca,'tickdir','out')
subplot(133),scatter(responseTot(:,2), responseTot(:,3), 24, 'filled','k');
hold on;
plot([0 1], [0 1], 'k--', 'LineWidth', 1.2); % y = x reference line
hold off;
xlabel('Mean Spike Response', 'FontWeight', 'bold');
ylabel('Mean Spike Response', 'FontWeight', 'bold');
title('Pull 2 vs Pull 3', 'FontWeight', 'bold');
set(gca, 'FontSize', 8, 'Box', 'off', 'GridAlpha', 0.4);
axis square;
box off,set(gca,'tickdir','out')
%%
bin_edges = 0:0.05:1; % Adjust bin size if desired
colors = [0 0 1; 0 0.5 0; 0.7 0 0]; % blue, green, red for pulls
temp = nan(length(selectivity_index),3);
figure; hold on;

for p = 1:num_pulls
    sel_for_pull = selectivity_index(best_pull == p); % Example: your selectivity indices for this pull
    temp(1:length(sel_for_pull),p) = sel_for_pull;
    % Bin and normalize
    counts = histcounts(sel_for_pull, bin_edges, 'Normalization', 'count');
    bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
    % Stairs plot
    stairs(bin_centers, counts, 'Color', colors(p,:), 'LineWidth', 2);
end

xlabel('Pull Preference (Selectivity) Index');
ylabel('Probability');
legend({'Pull 1', 'Pull 2', 'Pull 3'}, 'Location', 'Best');
title('Normalized Pull Preference Index Distribution (Stairs Plot)');
hold off;

%% FUNCTIONS
