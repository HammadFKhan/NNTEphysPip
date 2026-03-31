% This code analyzes the drift of neural state space initial conditions
% across trials in 3D. Using clustering based statistics and so on. It calculates the drift vector, magnitude, and
% direction for each trial transition. It then visualizes drift magnitude
% and directional angles, and plots the initial condition trajectory with
% drift vectors. This analysis helps quantify and interpret how neural
% trajectories change over time in all dimensions.
% Note that we concatenate trial conditions as to apply the same models for
% statistical comparison (ie. hit vs miss, hit vs FA, opto vs no opto)
if ~isfield(Spikes,'GPFA')
    Spikes = makeSpikeGPFA(Spikes);
    Spikes.GPFA.HitMiss.dat = [Spikes.GPFA.hit.dat,Spikes.GPFA.miss.dat];
    for n = 1:length(IntanBehaviour.hitTrace)%+1:length(Spikes.GPFA.BaselineOpto.dat) %fix trials
        Spikes.GPFA.HitMiss.dat(n).trialId = n;
    end
    Spikes.GPFA.MIHitFA.dat = [Spikes.GPFA.MIHit.dat,Spikes.GPFA.MIFA.dat];
    for n = length(IntanBehaviour.MIHitTrace)+1:length(Spikes.GPFA.MIHitFA.dat) %fix trials
        Spikes.GPFA.MIHitFA.dat(n).trialId = n;
    end
    Spikes.GPFA.HitEffort.dat = [Spikes.GPFA.MIHit.dat,Spikes.GPFA.effortperturb.dat];
    for n = 1:length(IntanBehaviour.MIHitTrace)%+1:length(Spikes.GPFA.BaselineOpto.dat) %fix trials
        Spikes.GPFA.HitEffort.dat(n).trialId = n;
    end
    
    addpath(genpath('C:\Users\khan332\Documents\GitHub\NeuralTraj'));
    addpath(genpath('mat_results'));
    if exist('mat_results','dir'),rmdir('mat_results','s'),end
    [Spikes.GPFA.resultHit,Spikes.GPFA.seqTrainHit] = gpfaAnalysis(Spikes.GPFA.hit.dat,1); %Run index
    [Spikes.GPFA.resultMiss,Spikes.GPFA.seqTrainMiss] = gpfaAnalysis(Spikes.GPFA.miss.dat,2); %Run index
    [Spikes.GPFA.resultMIHit,Spikes.GPFA.seqTrainMIHit] = gpfaAnalysis(Spikes.GPFA.MIHit.dat,3); %Run index
    % [Spikes.GPFA.resultMIFA,Spikes.GPFA.seqTrainMIFA] = gpfaAnalysis(Spikes.GPFA.MIFA.dat,4); %Run index
    [Spikes.GPFA.resultHitMiss,Spikes.GPFA.seqTrainHitMiss] = gpfaAnalysis(Spikes.GPFA.HitMiss.dat,5); %Run index
    % [Spikes.GPFA.resultMIHitFA,Spikes.GPFA.seqTrainMIHitFA] = gpfaAnalysis(Spikes.GPFA.MIHitFA.dat,6); %Run index
    [Spikes.GPFA.resultHitEffort,Spikes.GPFA.seqTrainHitEffort] = gpfaAnalysis(Spikes.GPFA.HitEffort.dat,7); %Run index
    close all
end
[neuralDynamics,waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);

%% Perform K-means clustering
% Load your data (replace `initCond` with actual variable if different)
% initCond = ... % n x 3 matrix
[c,allTrials] = sort_hit_effort(IntanBehaviour);
x = horzcat(squeeze(neuralDynamics.hiteffort.X(1,:,:)),squeeze(neuralDynamics.effort.X(1,:,:)));
y = horzcat(squeeze(neuralDynamics.hiteffort.X(2,:,:)),squeeze(neuralDynamics.effort.X(2,:,:)));
z = horzcat(squeeze(neuralDynamics.hiteffort.X(3,:,:)),squeeze(neuralDynamics.effort.X(3,:,:)));
% sort trials
x = x(45,c);
y = y(45,c);
z = z(45,c);
initCond = [x',y',z'];
% Choose the number of clusters, e.g., 2 or 3 (can be tuned or assessed later)
numClusters = 3;
[idx, C] = kmeans(initCond, numClusters, 'Replicates', 10);
% Assume idx is your cluster assignment vector (values 1~3 for three clusters)
cluster_colors = [1 0.3 0.3;    % red
                  0.8 .3 0.8;    % blue
                  0.0 0.45 1];   % purple

figure;
hold on
for k = 1:3
    scatter3(initCond(idx==k,1), ...
             initCond(idx==k,2), ...
             initCond(idx==k,3), ...
             50, cluster_colors(k,:), 'filled');
end

% Overlay cluster centers (assuming C is centers array)
scatter3(C(:,1), C(:,2), C(:,3), 100, 'k', 'x', 'LineWidth', 2);
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
grid on
hold off
view(30,30)
axis square
%%
t = 50;
[c,allTrials] = sort_hit_effort(IntanBehaviour);
x = horzcat(squeeze(neuralDynamics.hiteffort.X(1,:,:)),squeeze(neuralDynamics.effort.X(1,:,:)));
y = horzcat(squeeze(neuralDynamics.hiteffort.X(2,:,:)),squeeze(neuralDynamics.effort.X(2,:,:)));
z = horzcat(squeeze(neuralDynamics.hiteffort.X(3,:,:)),squeeze(neuralDynamics.effort.X(3,:,:)));
x = x(:,c);
y = y(:,c);
z = z(:,c);
figure,hold on
for trial = 1:size(x,2)
        % Plot a dot for the current time point based on effort
        if allTrials(2,trial)==1
            plot3(x(t,trial), y(t,trial), z(t,trial), 'o', 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', 'k');
        else
            plot3(x(t,trial), y(t,trial), z(t,trial), 'o', 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', 'k');
        end
end
view(30,30)
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
grid on
axis square
%%
figure
for k = 1:3
    x = horzcat(squeeze(neuralDynamics.hiteffort.X(1,:,:)),squeeze(neuralDynamics.effort.X(1,:,:)));
    y = horzcat(squeeze(neuralDynamics.hiteffort.X(2,:,:)),squeeze(neuralDynamics.effort.X(2,:,:)));
    z = horzcat(squeeze(neuralDynamics.hiteffort.X(3,:,:)),squeeze(neuralDynamics.effort.X(3,:,:)));
    x = mean(x(:,idx==k),2);
    y = mean(y(:,idx==k),2);
    z = mean(z(:,idx==k),2);
    plot3(x,y,z),hold on
end
%%
[c,allTrials] = sort_hit_effort(IntanBehaviour);
x = horzcat(squeeze(neuralDynamics.hiteffort.X(1,:,:)),squeeze(neuralDynamics.effort.X(1,:,:)));
y = horzcat(squeeze(neuralDynamics.hiteffort.X(2,:,:)),squeeze(neuralDynamics.effort.X(2,:,:)));
z = horzcat(squeeze(neuralDynamics.hiteffort.X(3,:,:)),squeeze(neuralDynamics.effort.X(3,:,:)));
% sort trials
x = x(:,c);
y = y(:,c);
z = z(:,c);
timeEnd = 250;
nTrials = size(x, 2);

v = VideoWriter('D:\SQLever\neural_conditionsDLSDay17renew2.avi'); % Name your output file
v.FrameRate = 10; % Set the frame rate
open(v);
initial_azimuth = 30;
elevation = 45;

figure('Color', 'w');
hold on
axis tight
view(initial_azimuth,elevation)
xlabel('X')
ylabel('Y')
zlabel('Z')
hitTrials = size(neuralDynamics.hiteffort.X,3);
for t = 50
    clf; % Clear the figure each frame
    hold on
    % Plot each trajectory up to time t
    for trial = 1:size(x,2)
        plot3(x(t,1:trial), y(t,1:trial), z(t,1:trial), 'Color', [0 0 0 0.4]);
        % Plot a dot for the current time point based on effort
        if allTrials(2,trial)==1
            plot3(x(t,trial), y(t,trial), z(t,trial), 'o', 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', 'k');
        else
            plot3(x(t,trial), y(t,trial), z(t,trial), 'o', 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', 'k');
        end

        title(['Initial conditions on trial = ', num2str(trial)]);
        axis([min(x(:)), max(x(:)), min(y(:)), max(y(:)), min(z(:)), max(z(:))]);
        % Calculate current azimuth angle for rotation
        current_azimuth = mod(initial_azimuth + 0.5*trial, 360);
        view(current_azimuth, elevation);
        grid on
        drawnow
        % Capture the frame and write to video
        frame = getframe(gcf);
        writeVideo(v, frame);
    end
end

close(v);
%% Proportion of effort and non effort trials in each cluster
isnoeffort = allTrials(2,:);
% sort sorted cluser based on chrnology so we can map onto the allTrials
% index
sortedclusterId = idx;
for k = 1:3
% return allTrial value of zero or 1
clusterTrials = isnoeffort(sortedclusterId==k); % returns trial number of that cluster
% sum of the cluster Trial is the number of non effort trials
effortproportion(k) = 1-(sum(clusterTrials)/length(clusterTrials)); % qik mathss
end
% plot out a bar plot
figure,bar(effortproportion','stacked')
%% Calculate Mahalanobis distance of clusters
n = size(initCond, 1);
mahalDists = zeros(n, numClusters);
for k = 1:numClusters
    mu = C(k, :);
    sigma = cov(initCond(idx==k, :));      % Covariance of cluster k
    for i = 1:n
        mahalDists(i, k) = sqrt((initCond(i,:) - mu) / sigma * (initCond(i,:) - mu)'); % Mahalanobis distance formula
    end
end
% Visualize or output as needed
[sortedDist, trialOrder] = sort(mahalDists(:,1), 'descend');
mahalDists = mahalDists(trialOrder, :); % Now rows are ordered by their dist to Cluster 1
nTrials = size(mahalDists,1);


addpath(genpath('C:\Users\khan332\Documents\GitHub\slanCM'));
figure,hold on
plotNiceBars(mahalDists)
ylim([0 ceil(max(mahalDists,[],'all'))])
xlabel('Cluster');
ylabel('Mahalanobis Distance from Center');
title('Within-cluster Mahalanobis Distances');
normDist = (sortedDist - min(sortedDist)) / (max(sortedDist) - min(sortedDist));
cmap = slanCM('RdBu',nTrials); % or any other colormap
trialColors = cmap(round(normDist * (size(cmap,1)-1))+1, :);
colormap(slanCM('RdBu'))
c = colorbar;
c.Label.String = 'Sorted Mahalanobis Distance';
c.Label.FontSize = 8;
% Optionally set ticks to match the real value range:
c.Ticks = [0 0.5 1];
c.TickLabels = {num2str(min(sortedDist)), num2str(mean(sortedDist)), num2str(max(sortedDist))};
axis square
%% Silohette calculation
n = size(initCond,1);
pairwiseMahal = zeros(n,n);

% Use overall covariance for simplicity
Sigma = cov(initCond);

for i = 1:n
    for j = 1:n
        diffs = initCond(i,:) - initCond(j,:);
        pairwiseMahal(i,j) = sqrt(diffs / Sigma * diffs');
    end
end
silo = zeros(n,1);

for i = 1:n
    myCluster = idx(i);
    sameInds = find(idx == myCluster & (1:n)' ~= i);
    otherClusters = setdiff(unique(idx), myCluster);
    
    % a(i): mean Mahalanobis distance to same cluster
    a_i = mean(pairwiseMahal(i, sameInds));
    
    % b(i): minimum mean Mahalanobis distance to other clusters
    b_i = inf;
    for k = otherClusters'
        kInds = find(idx == k);
        b_ik = mean(pairwiseMahal(i, kInds));
        if b_ik < b_i
            b_i = b_ik;
        end
    end
    
    silo(i) = (b_i - a_i) / max(a_i, b_i);
end

s = silhouette(initCond,idx);

% Plot it
figure;
h = histogram(silo, 11,'Normalization', 'probability', 'FaceColor', [0.5 0.8 1], 'EdgeColor','none');
set(gca, ...
    'TickDir', 'out', ...
    'Box', 'off', ...
    'FontSize', 14, ...
    'LineWidth', 1.5);
axis square;
xlabel('Silhouette Value', 'FontSize', 16);
ylabel('Probability', 'FontSize', 16);
title('Silhouette', 'FontSize', 16, 'FontWeight', 'normal');
% Set consistent limits for clarity
xlim([-0.5, 0.75]); % Adjust as needed for your data
ylim([0, max(h.Values)*1.1]);
% Compute and plot mean
m = median(silo);
skewness(silo)
yl = ylim;
hold on;
xline(m, '--k', ['Median = ' num2str(m, '%.2f')], ...
    'LineWidth', 2, ...
    'LabelOrientation', 'horizontal', ...
    'LabelHorizontalAlignment', 'center', ...
    'LabelVerticalAlignment', 'top', ...
    'FontSize', 8, ...
    'Color', [0.3 0.3 0.3]);
hold off;

%% Drift metrics

% Assume initCond is N_trials x dim
N = size(initCond, 1);
driftVec = [diff(initCond)]; % (N-1) x dim, each row is drift vector to next trial
driftMag = sqrt(sum(driftVec.^2, 2)); % Magnitude of drift for each transition
dirCosine = driftVec ./ driftMag; % (N-1) x 3, each row gives [cos_alpha, cos_beta, cos_gamma] for that drift
% XY plane angle
theta_xy = atan2(driftVec(:,2), driftVec(:,1)); % angle in XY plane
% YZ plane angle
theta_yz = atan2(driftVec(:,3), driftVec(:,2)); % angle in YZ plane
% XZ plane angle
theta_xz = atan2(driftVec(:,3), driftVec(:,1)); % angle in XZ plane
theta_xy_deg = rad2deg(theta_xy);
theta_yz_deg = rad2deg(theta_yz);
theta_xz_deg = rad2deg(theta_xz);

% Spherical coordinates: azimuth and elevation
dx = driftVec(:,1);
dy = driftVec(:,2);
dz = driftVec(:,3);

azimuth = atan2(dy, dx); % Angle in XY plane
elevation = atan2(dz, sqrt(dx.^2 + dy.^2)); % Angle from XY plane

% Convert to degrees for easier interpretation
azimuth_deg = rad2deg(azimuth);
elevation_deg = rad2deg(elevation);
%%
figure;
% Top subplot: Drift magnitude (smoothed)
subplot(2,1,1);
plot(smoothdata(driftMag,'movmean',7), 'k.','MarkerSize', 15);
ylabel('Drift Magnitude');
xlabel('Trial Transition');
box off; % Remove top/right borders for clean look
set(gca, 'TickDir', 'out');
ylim([0 1])
% Bottom subplot: Overlay azimuth and elevation
subplot(2,1,2);
plot(smoothdata(azimuth_deg,'movmean',5),'.', 'MarkerSize', 15,'Color', [0.25 0.55 0.88]); hold on;
plot(smoothdata(elevation_deg,'movmean',5),'.', 'MarkerSize', 15,'Color', [0.65 0.2 0.55]); % 'm' for magenta (similar to your example)
legend({'Azimuth phase', 'Elevation phase'}, 'Location', 'best');
ylabel('Phase (deg)');
xlabel('Trial Transition');
box off; set(gca, 'TickDir', 'out');

sgtitle('Trajectory Drift and Angular Phase Response');
%%
% Assuming driftMag (length N-1) and silhouette_scores (length N)
% Align vectors (skip first silhouette for direct correspondence)
X = zscore(driftMag); 
Y = zscore(s);

% Run linear regression and display statistics
mdl = fitlm(X,Y); % Regression model
disp(mdl);

% Plot
figure;
scatter(X, Y, 40, 'filled', 'MarkerFaceColor', [0.25 0.55 0.88]); % Blue, match prior visual style
hold on;
plot(X, mdl.Fitted, '-r', 'LineWidth', 2); % Regression line in red
xlabel('Drift Magnitude');
ylabel('Silhouette Score');
title('Drift Magnitude vs Cluster Membership');
grid on;
legend({'Data', 'Regression fit'}, 'Location', 'best');
%% Cluster label spikes
spike_thresh = prctile(driftMag, 90);
spike_trials = find(driftMag > spike_thresh); % Index of trials where drift is high
is_switch = idx(spike_trials) ~= idx(spike_trials+1);


figure; hold on;
plot(driftMag, 'k-', 'LineWidth', 1.5);
scatter(spike_trials(~is_switch), driftMag(spike_trials(~is_switch)), 60, 'b', 'filled'); % High drift, same cluster
scatter(spike_trials(is_switch), driftMag(spike_trials(is_switch)), 60, 'r', 'filled'); % High drift, cluster switch
ylabel('Drift Magnitude');
xlabel('Trial');
legend({'Drift', 'High Drift, same cluster', 'High Drift, switched cluster'}, 'Location', 'best');
title('Drift Magnitude and Cluster Switching Events');

cluster_colors = [1 0.3 0.3;    % red
                  0.8 .3 0.8;    % blue
                  0.0 0.45 1];   % purple

% Create an N x 3 color matrix for each point
N = numel(idx);
point_colors = cluster_colors(idx, :);

figure; hold on;
scatter(1:N, idx, 40, point_colors, 'filled');

for i = 1:numel(spike_trials)
    plot([spike_trials(i) spike_trials(i)], ylim, '--', 'Color', [0.8 0.8 0.8]); % Vertical drift spike
end
ylabel('Cluster Label');
xlabel('Trial');
title('Cluster Labels with Drift Spike Events');
%% Plot out significant initial condition jumps as a function of effort trials
% Assume initCond is N_trials x dim
alignedPerturbations = allTrials(2,:)==1;
initCondCatch = initCond(alignedPerturbations,:);

idxCatch = idx(alignedPerturbations);

N = size(initCondCatch, 1);
driftVec = [diff(initCondCatch)]; % (N-1) x dim, each row is drift vector to next trial
driftMag = sqrt(sum(driftVec.^2, 2)); % Magnitude of drift for each transition
dirCosine = driftVec ./ driftMag; % (N-1) x 3, each row gives [cos_alpha, cos_beta, cos_gamma] for that drift
% XY plane angle
theta_xy = atan2(driftVec(:,2), driftVec(:,1)); % angle in XY plane
% YZ plane angle
theta_yz = atan2(driftVec(:,3), driftVec(:,2)); % angle in YZ plane
% XZ plane angle
theta_xz = atan2(driftVec(:,3), driftVec(:,1)); % angle in XZ plane
theta_xy_deg = rad2deg(theta_xy);
theta_yz_deg = rad2deg(theta_yz);
theta_xz_deg = rad2deg(theta_xz);

% Spherical coordinates: azimuth and elevation
dx = driftVec(:,1);
dy = driftVec(:,2);
dz = driftVec(:,3);

azimuth = atan2(dy, dx); % Angle in XY plane
elevation = atan2(dz, sqrt(dx.^2 + dy.^2)); % Angle from XY plane

% Convert to degrees for easier interpretation
azimuth_deg = rad2deg(azimuth);
elevation_deg = rad2deg(elevation);

driftMagNormative = driftMag;
silhouetteNormative = s(2:end);
% assert(length(silhouetteNormative)==size(M1neuralDynamics.hiteffort.X,3))

spike_thresh = prctile(driftMagNormative, 0); % Use 0 for all transitions, adjust if you want quantiles
significant_jumps = find(driftMagNormative > spike_thresh);
same_cluster = diff(idxCatch)==0;
% same_cluster = idx(significant_jumps) == idx(significant_jumps + 1);
switch_cluster = ~same_cluster;

% Format for plotNiceBars: 2 columns, pad with NaNs for unequal group sizes
maxN = max(sum(same_cluster), sum(switch_cluster));
totSilhouette = nan(maxN,2);
totDrift = nan(maxN,2);

% Fill columns
totSilhouette(1:sum(same_cluster),1) = silhouetteNormative(same_cluster);
totSilhouette(1:sum(switch_cluster),2) = silhouetteNormative(switch_cluster);

totDrift(1:sum(same_cluster),1) = driftMagNormative(same_cluster);
totDrift(1:sum(switch_cluster),2) = driftMagNormative(switch_cluster);

% Proportion data (for bar plot, not for plotNiceBars)
props = [sum(same_cluster), sum(switch_cluster)] / numel(significant_jumps);

% --- Plotting ---
figure;
bar(categorical({'Stayed','Switched'}), props, 'FaceColor',[0.8 0.8 0.8],'EdgeColor','k');
ylabel('Proportion of Transitions');
title('Proportion: Stay vs Switch Cluster');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 14, 'LineWidth', 1.5);

figure;hold on
plotNiceBars(totSilhouette);
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 14, 'LineWidth', 1.5);
title('Silhouette Score by Cluster Transition');
set(gca, 'XTickLabel', {'Stayed','Switched'});
ylabel('Silhouette Score');
ylim([0 1])
figure;hold on
plotNiceBars(totDrift);
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 14, 'LineWidth', 1.5);
title('Drift Magnitude by Cluster Transition');
set(gca, 'XTickLabel', {'Stayed','Switched'});
ylabel('Drift Magnitude');


% Elevation and azimuth deg
nBins = floor(length(azimuth_deg)/10);
edges = linspace(0, 2*pi, nBins+1); % Bin edges for full circle

% Assume azimuth_rad and elevation_rad are already in radians,
% and driftMag is the same length (one per trial)
azimuth_rad = mod(deg2rad(azimuth_deg), 2*pi);
elevation_rad = mod(deg2rad(elevation_deg), 2*pi);

% Bin centers
bin_centers = edges(1:end-1) + diff(edges)/2;

% Calculate mean drift magnitude per angle bin (azimuth)
az_drift_mean = zeros(1, nBins);
for i = 1:nBins
    in_bin = azimuth_rad >= edges(i) & azimuth_rad < edges(i+1);
    if sum(in_bin)>0
        az_drift_mean(i) = mean(driftMag(in_bin));
    end
end

% Same for elevation
el_drift_mean = zeros(1, nBins);
for i = 1:nBins
    in_bin = elevation_rad >= edges(i) & elevation_rad < edges(i+1);
    if sum(in_bin)>0
        el_drift_mean(i) = mean(driftMag(in_bin));
    end
end

% Interpolate the binned line data for smoother curves
nInterp = 150; % Number of points for smooth line
theta_interp = linspace(0, 2*pi, nInterp);

az_drift_smooth = smoothdata(interp1(bin_centers, az_drift_mean, theta_interp, 'linear'),'gaussian',5);
el_drift_smooth = smoothdata(interp1(bin_centers, el_drift_mean, theta_interp, 'linear'),'gaussian',5);

figure;
pax = polaraxes;
hold on;

% Raw (unbinned) scatter for azimuth
% polarscatter(azimuth_rad, smoothdata(driftMag,'gaussian',1), 18, 'k', 'filled');
% Raw (unbinned) scatter for elevation (optional)
% polarscatter(elevation_rad, driftMag, 18, [0.27 0.69 0.98], 'filled');

% Smoothed lines for binned means
polarplot(theta_interp, az_drift_smooth, 'Color', [0.25 0.55 0.88], 'LineWidth', 2);
polarplot(theta_interp, el_drift_smooth, 'Color', [0.65 0.2 0.55], 'LineWidth', 2);

% Format
set(gca, 'FontSize', 12);
legend({'Azimuth Drift Mean','Elevation Drift Mean'}, 'Location', 'northeastoutside');
hold off;

%% Now plot out as a function of break trials
alignedPerturbations = allTrials(2,:)==0;
initCondCatch = initCond(alignedPerturbations,:);

idxCatch = idx(alignedPerturbations);

N = size(initCondCatch, 1);
driftVec = [diff(initCondCatch)]; % (N-1) x dim, each row is drift vector to next trial
driftMag = sqrt(sum(driftVec.^2, 2)); % Magnitude of drift for each transition
dirCosine = driftVec ./ driftMag; % (N-1) x 3, each row gives [cos_alpha, cos_beta, cos_gamma] for that drift
% XY plane angle
theta_xy = atan2(driftVec(:,2), driftVec(:,1)); % angle in XY plane
% YZ plane angle
theta_yz = atan2(driftVec(:,3), driftVec(:,2)); % angle in YZ plane
% XZ plane angle
theta_xz = atan2(driftVec(:,3), driftVec(:,1)); % angle in XZ plane
theta_xy_deg = rad2deg(theta_xy);
theta_yz_deg = rad2deg(theta_yz);
theta_xz_deg = rad2deg(theta_xz);

% Spherical coordinates: azimuth and elevation
dx = driftVec(:,1);
dy = driftVec(:,2);
dz = driftVec(:,3);

azimuth = atan2(dy, dx); % Angle in XY plane
elevation = atan2(dz, sqrt(dx.^2 + dy.^2)); % Angle from XY plane

% Convert to degrees for easier interpretation
azimuth_deg = rad2deg(azimuth);
elevation_deg = rad2deg(elevation);

driftMagNormative = driftMag;
silhouetteNormative = s(2:end);
% assert(length(silhouetteNormative)==size(M1neuralDynamics.hiteffort.X,3))

spike_thresh = prctile(driftMagNormative, 0); % Use 0 for all transitions, adjust if you want quantiles
significant_jumps = find(driftMagNormative > spike_thresh);
same_cluster = diff(idxCatch)==0;
% same_cluster = idx(significant_jumps) == idx(significant_jumps + 1);
switch_cluster = ~same_cluster;

% Format for plotNiceBars: 2 columns, pad with NaNs for unequal group sizes
maxN = max(sum(same_cluster), sum(switch_cluster));
totSilhouette = nan(maxN,2);
totDrift = nan(maxN,2);

% Fill columns
totSilhouette(1:sum(same_cluster),1) = silhouetteNormative(same_cluster);
totSilhouette(1:sum(switch_cluster),2) = silhouetteNormative(switch_cluster);

totDrift(1:sum(same_cluster),1) = driftMagNormative(same_cluster);
totDrift(1:sum(switch_cluster),2) = driftMagNormative(switch_cluster);

% Proportion data (for bar plot, not for plotNiceBars)
props = [sum(same_cluster), sum(switch_cluster)] / numel(significant_jumps);

% --- Plotting ---
figure;
bar(categorical({'Stayed','Switched'}), props, 'FaceColor',[0.8 0.8 0.8],'EdgeColor','k');
ylabel('Proportion of Transitions');
title('Proportion: Stay vs Switch Cluster');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 14, 'LineWidth', 1.5);

figure;hold on
plotNiceBars(totSilhouette);
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 14, 'LineWidth', 1.5);
title('Silhouette Score by Cluster Transition');
set(gca, 'XTickLabel', {'Stayed','Switched'});
ylabel('Silhouette Score');

figure;hold on
plotNiceBars(totDrift);
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 14, 'LineWidth', 1.5);
title('Drift Magnitude by Cluster Transition');
set(gca, 'XTickLabel', {'Stayed','Switched'});
ylabel('Drift Magnitude');
 
% Elevation and azimuth deg
nBins = floor(length(azimuth_deg)/10);
edges = linspace(0, 2*pi, nBins+1); % Bin edges for full circle

% Assume azimuth_rad and elevation_rad are already in radians,
% and driftMag is the same length (one per trial)
azimuth_rad = mod(deg2rad(azimuth_deg), 2*pi);
elevation_rad = mod(deg2rad(elevation_deg), 2*pi);

% Bin centers
bin_centers = edges(1:end-1) + diff(edges)/2;

% Calculate mean drift magnitude per angle bin (azimuth)
az_drift_mean = zeros(1, nBins);
for i = 1:nBins
    in_bin = azimuth_rad >= edges(i) & azimuth_rad < edges(i+1);
    if sum(in_bin)>0
        az_drift_mean(i) = mean(driftMag(in_bin));
    end
end

% Same for elevation
el_drift_mean = zeros(1, nBins);
for i = 1:nBins
    in_bin = elevation_rad >= edges(i) & elevation_rad < edges(i+1);
    if sum(in_bin)>0
        el_drift_mean(i) = mean(driftMag(in_bin));
    end
end

% Interpolate the binned line data for smoother curves
nInterp = 150; % Number of points for smooth line
theta_interp = linspace(0, 2*pi, nInterp);

az_drift_smooth = smoothdata(interp1(bin_centers, az_drift_mean, theta_interp, 'linear'),'gaussian',5);
el_drift_smooth = smoothdata(interp1(bin_centers, el_drift_mean, theta_interp, 'linear'),'gaussian',5);

figure;
pax = polaraxes;
hold on;

% Raw (unbinned) scatter for azimuth
% polarscatter(azimuth_rad, smoothdata(driftMag,'gaussian',1), 18, 'k', 'filled');
% Raw (unbinned) scatter for elevation (optional)
% polarscatter(elevation_rad, driftMag, 18, [0.27 0.69 0.98], 'filled');

% Smoothed lines for binned means
polarplot(theta_interp, az_drift_smooth, 'Color', [0.25 0.55 0.88], 'LineWidth', 2);
polarplot(theta_interp, el_drift_smooth, 'Color', [0.65 0.2 0.55], 'LineWidth', 2);

% Format
set(gca, 'FontSize', 12);
legend({'Azimuth Drift Mean','Elevation Drift Mean'}, 'Location', 'northeastoutside');
hold off;
%%
function plotNiceBars(totData)
means = nanmean(totData);          % Bar heights
sems = nanstd(totData) ./ sqrt(size(totData,1));   % Error bar (standard error)
b = bar(means, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k'); % Gray bars with black edge

% Overlay error bars
% errorbar(1:size(totData,2), means, sems, 'k', 'LineStyle', 'none', 'LineWidth', 1);

% Overlay individual jittered points
xjitter = randn(size(totData))*0.03; % Controls point jitter
for i = 1:size(totData,2)
    scatter(i + xjitter(:,min(i,2)), totData(:,i), 50, '.',...
        'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4);
end

nTrials = size(totData,1);
% cm = slanCM('RdBu',nTrials); % Or use your favorite colormap
% % % Draw paired lines between columns 1 and 2
% for j = 1:size(totData,1)
%     if size(totData,2) >= 3  % If there are at least 3 columns
%         xvals = [1 + xjitter(j,1), 2 + xjitter(j,2), 3 + xjitter(j,3)];
%         yvals = [totData(j,1),    totData(j,2),    totData(j,3)];
%         plot(xvals, yvals, '-', 'Color', cm(j,:), 'LineWidth', 1);
%     else % Connect just columns 1 and 2
%         xvals = [1 + xjitter(j,1), 2 + xjitter(j,2)];
%         yvals = [totData(j,1),    totData(j,2)];
%         plot(xvals, yvals, '-', 'Color', cm(j,:), 'LineWidth', 1);
%     end
% end

% Style similar to image
set(gca, 'XTick', 1:size(totData,2), 'XTickLabel', {'Second Pull', 'Third Pull', 'Polymer', 'Late'}, ...
    'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
ylabel('IPI (s)');
% ylim([0 2]);

hold off;

%%% RUN STATS
[p, tbl, stats] = anova1(totData, [], 'off'); % columns as groups
results = multcompare(stats, 'Display', 'off'); % Pairwise comparisons

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
