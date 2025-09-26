% Warp response pooled together for M1 and DLS
M1Sq= struct();
DLSSq = struct();
files = dir(fullfile('D:\SequenceProject\WarpedSpikes\M1\','*.mat'));
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    [warpedSpks.GPFA.p1neuralDynamics] = getGPFASq(warpedSpks.pull1A.warpSpikes,IntanBehaviour);
    [warpedSpks.GPFA.p2neuralDynamics] = getGPFASq(warpedSpks.pull2A.warpSpikes,IntanBehaviour);
    [warpedSpks.GPFA.p3neuralDynamics] = getGPFASq(warpedSpks.pull3A.warpSpikes,IntanBehaviour);
    M1Sq(fileNum).warpedSpikes = warpedSpks;
    M1Sq(fileNum).IntanBehaviour = IntanBehaviour;
    M1Sq(fileNum).filename = files(fileNum).name;
    prepRSLDS(warpedSpks.warpSpikes,fullfile(files(fileNum).folder,files(fileNum).name));
end

files = dir(fullfile('D:\SequenceProject\WarpedSpikes\DLS\','*.mat'));
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    [warpedSpks.GPFA.p1neuralDynamics] = getGPFASq(warpedSpks.pull1A.warpSpikes,IntanBehaviour);
    [warpedSpks.GPFA.p2neuralDynamics] = getGPFASq(warpedSpks.pull2A.warpSpikes,IntanBehaviour);
    [warpedSpks.GPFA.p3neuralDynamics] = getGPFASq(warpedSpks.pull3A.warpSpikes,IntanBehaviour);
    DLSSq(fileNum).warpedSpikes = warpedSpks;
    DLSSq(fileNum).IntanBehaviour = IntanBehaviour;
    DLSSq(fileNum).filename = files(fileNum).name;
    prepRSLDS(warpedSpks.warpSpikes,fullfile(files(fileNum).folder,files(fileNum).name));
end
%% Representative M1 trace
neuralDynamics = M1Sq(3).warpedSpikes.GPFA.p3neuralDynamics;
time = M1Sq(3).warpedSpikes.warpedTime(:,3);
time = linspace(time(1),time(end),size(neuralDynamics.hitOnly.X,2));
figure,
subplot(311),plot(time,squeeze(neuralDynamics.hitOnly.X(1,:,:)),'Color',[0 0 0 0.4]),hold on,box off,set(gca,'fontsize',18)
subplot(312),plot(time,squeeze(neuralDynamics.hitOnly.X(2,:,:)),'Color',[0 0 0 0.4]),hold on,box off,set(gca,'fontsize',18)
subplot(313),plot(time,squeeze(neuralDynamics.hitOnly.X(3,:,:)),'Color',[0 0 0 0.4]),hold on,box off,set(gca,'fontsize',18)
pullIndex = vertcat(IntanBehaviour.hitTrace.pullCount);
x = squeeze(neuralDynamics.hitOnly.X(1,:,:));
y = squeeze(neuralDynamics.hitOnly.X(2,:,:));
z = squeeze(neuralDynamics.hitOnly.X(3,:,:));
% figure,hold on
% for n = 1:size(neuralDynamics.hitOnly.X,3)
% plot3(x(:,n),y(:,n),z(:,n),'color',[0 0 0 0.4])
% end
colors = [12,188,187;183,13,180]/255;
timeIndex = linspace(1,size(x,1),5001);
xp = mean(x,2);
yp = mean(y,2);
zp = mean(z,2);
figure
plot3(xp,yp,zp,'color',colors(2,:)),hold on
pIm = floor(mean(pullIndex));
pIm = [1 pIm 3499]; %add reward 
for n = 1:length(pIm)
    id = floor(timeIndex(pIm(n)));
    scatter3(xp(id),yp(id),zp(id),20,'filled'),hold on
end
%% Representative DLS trace
neuralDynamics = DLSSq(5).warpedSpikes.GPFA.p3neuralDynamics;
time = DLSSq(5).warpedSpikes.warpedTime(:,3);
time = linspace(time(1),time(end),size(neuralDynamics.hitOnly.X,2));
figure,
subplot(311),plot(time,squeeze(neuralDynamics.hitOnly.X(1,:,:)),'Color',[0 0 0 0.4]),hold on,box off,set(gca,'fontsize',18)
subplot(312),plot(time,squeeze(neuralDynamics.hitOnly.X(2,:,:)),'Color',[0 0 0 0.4]),hold on,box off,set(gca,'fontsize',18)
subplot(313),plot(time,squeeze(neuralDynamics.hitOnly.X(3,:,:)),'Color',[0 0 0 0.4]),hold on,box off,set(gca,'fontsize',18)
pullIndex = vertcat(IntanBehaviour.hitTrace.pullCount);
x = squeeze(neuralDynamics.hitOnly.X(1,:,:));
y = squeeze(neuralDynamics.hitOnly.X(2,:,:));
z = squeeze(neuralDynamics.hitOnly.X(3,:,:));
% figure,hold on
% for n = 1:size(neuralDynamics.hitOnly.X,3)
% plot3(x(:,n),y(:,n),z(:,n),'color',[0 0 0 0.4])
% end
colors = [12,188,187;183,13,180]/255;
timeIndex = linspace(1,size(x,1),5001);
xp = mean(x,2);
yp = mean(y,2);
zp = mean(z,2);
figure
plot3(xp,yp,zp,'color',colors(1,:)),hold on
pIm = floor(mean(pullIndex));
pIm = [1 pIm 3499]; %add reward 
for n = 1:length(pIm)
    id = floor(timeIndex(pIm(n)));
    scatter3(xp(id),yp(id),zp(id),20,'filled'),hold on
end
%%
% figure,hold on
% for n = 1:200
%     plot(time,squeeze(M1neuralDynamics.hit.speed.speed(1,:,n)),'color',[0 0 0 0.8])
% %     for p = 1:length(IntanBehaviour.hitTrace(n).pullCount)
% %         xline((IntanBehaviour.hitTrace(n).pullCount(p)-3500)/1000,'r')
% %     end
% end
neuralDynamics = DLSSq(5).warpedSpikes.GPFA.p3neuralDynamics;
time = linspace(time(1),time(end),size(neuralDynamics.hitOnly.X,2));
colors = [12,188,187;183,13,180]/255;
speedTotBaseline = smoothdata(squeeze(neuralDynamics.hitOnly.speed.speed(1,2:end,:)),1,'movmean',10);
figure,hold on
plot(time(2:end),mean(speedTotBaseline,2),'color',colors(1,:),'linewidth',2),hold on
plot(time(2:end),mean(speedTotBaseline,2)+std(speedTotBaseline,[],2)/(sqrt(size(speedTotBaseline,2))),'color',colors(1,:),'linewidth',2)
plot(time(2:end),mean(speedTotBaseline,2)-std(speedTotBaseline,[],2)/(sqrt(size(speedTotBaseline,2))),'color',colors(1,:),'linewidth',2)
set(gca,'tickdir','out'),box off, axis square
xlim([-3.5,1.5])

neuralDynamics = M1Sq(3).warpedSpikes.GPFA.p1neuralDynamics;
time = linspace(time(1),time(end),size(neuralDynamics.hitOnly.X,2));
colors = [12,188,187;183,13,180]/255;
speedTotBaseline = smoothdata(squeeze(neuralDynamics.hitOnly.speed.speed(1,2:end,:)),1,'gaussian',5);
figure,hold on
plot(time(2:end),mean(speedTotBaseline,2),'color',colors(2,:),'linewidth',2),hold on
plot(time(2:end),mean(speedTotBaseline,2)+std(speedTotBaseline,[],2)/(sqrt(size(speedTotBaseline,2))),'color',colors(2,:),'linewidth',2)
plot(time(2:end),mean(speedTotBaseline,2)-std(speedTotBaseline,[],2)/(sqrt(size(speedTotBaseline,2))),'color',colors(2,:),'linewidth',2)
set(gca,'tickdir','out'),box off, axis square
xlim([-3.5,1.5])
%% Lets calculate the average speed modulation based on aligned spikes
% Speed modulation can be the modulation based on the change in baseline. 
% Here we calculate the maximum speed at the time of pull from baseline
% (-200 ms) and the pull response (+100 ms);
M1modulation = cell(1,length(M1Sq));
for n = 1:length(M1Sq)
    neuralDynamics = M1Sq(n);
    M1modulation{n} = getModulation(neuralDynamics);
end
%%
% M1modulation = vertcat(M1modulation{:});
% figure,
% customBarplot(M1modulation)
% ylabel('Speed Modulation'),axis square
% set(gca,'tickdir','out'),box off

%M1modulation = [M1modulation,cellfun(@(x) x(1:floor(length(x))/2,:),M1modulation,'UniformOutput',false)];
cellCount = numel(M1modulation);                % Number of cells
pulls = size(M1modulation{1}, 2);               % Number of pulls

pullMeans = nan(cellCount, pulls);              % Preallocate for cell means

for i = 1:cellCount
    pullMeans(i, :) = mean(M1modulation{i}, 1); % Mean across trials
end

mean_vals = mean(pullMeans, 1);
sem_vals = std(pullMeans, 0, 1) / sqrt(cellCount);

figure;
hold on;
colors = lines(cellCount);
% Plot individual cell lines and dots
for i = 1:cellCount
    plot(1:pulls, pullMeans(i, :), '-o', ...
        'Color', colors(i, :), ... 
        'MarkerFaceColor', colors(i, :), ...
        'MarkerEdgeColor', colors(i, :), ...
        'LineWidth', 1.5, 'MarkerSize', 7);
end
% Overlay the barplot (mean ± SEM)
bar_handle = bar(1:pulls, mean_vals, 0.7, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
for b = 1:pulls
    bar_handle.FaceColor = [0.4 0.4 0.4]; % Optional: neutral bar color
end
xlabel('Pull');
ylabel('Average Modulation');
set(gca, 'XTick', 1:pulls,'tickdir','out');
hold off;
xlim([0 4])
axis square
%%
DLSmodulation = cell(1,length(DLSSq));
for n = 1:length(DLSSq)
    neuralDynamics = DLSSq(n);
    DLSmodulation{n} = getModulation(neuralDynamics);
end
%%
cellCount = numel(DLSmodulation);                % Number of cells
pulls = size(DLSmodulation{1}, 2);               % Number of pulls

pullMeans = nan(cellCount, pulls);              % Preallocate for cell means

for i = 1:cellCount
    pullMeans(i, :) = mean(DLSmodulation{i}, 1); % Mean across trials
end

mean_vals = mean(pullMeans, 1);
sem_vals = std(pullMeans, 0, 1) / sqrt(cellCount);

figure;
hold on;
colors = lines(cellCount);
% Plot individual cell lines and dots
for i = 1:cellCount
    plot(1:pulls, pullMeans(i, :), '-o', ...
        'Color', colors(i, :), ... 
        'MarkerFaceColor', colors(i, :), ...
        'MarkerEdgeColor', colors(i, :), ...
        'LineWidth', 1.5, 'MarkerSize', 7);
end
% Overlay the barplot (mean ± SEM)
bar_handle = bar(1:pulls, mean_vals, 0.7, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
for b = 1:pulls
    bar_handle.FaceColor = [0.4 0.4 0.4]; % Optional: neutral bar color
end
xlabel('Pull');
ylabel('Average Modulation');
set(gca, 'XTick', 1:pulls,'tickdir','out');
hold off;
xlim([0 4])
axis square
%% Compute tangle
QM = [];
QD = [];
for n = 1:length(M1Sq)
Qm1 = computeTangling(M1Sq(n).warpedSpikes.GPFA.p1neuralDynamics.hitOnly.X,0.0001);
Qm1 = max(Qm1,[],1);
QM{n}(:,1) = Qm1;

Qm1 = computeTangling(M1Sq(n).warpedSpikes.GPFA.p2neuralDynamics.hitOnly.X,0.0001);
Qm1 = max(Qm1,[],1);
QM{n}(:,2) = Qm1;

Qm1 = computeTangling(M1Sq(n).warpedSpikes.GPFA.p3neuralDynamics.hitOnly.X,0.0001);
Qm1 = max(Qm1,[],1);
QM{n}(:,3) = Qm1;
end

for n = 1:length(DLSSq)
Qdls = computeTangling(DLSSq(n).warpedSpikes.GPFA.p1neuralDynamics.hitOnly.X,0.0001);
Qdls = max(Qdls,[],1);
Qm1 = max(Qm1,[],1);
QD{n}(:,1) = Qdls;

Qdls = computeTangling(DLSSq(n).warpedSpikes.GPFA.p2neuralDynamics.hitOnly.X,0.0001);
Qdls = max(Qdls,[],1);
QD{n}(:,2) = Qdls;

Qdls = computeTangling(DLSSq(n).warpedSpikes.GPFA.p3neuralDynamics.hitOnly.X,0.0001);
Qdls = max(Qdls,[],1);
QD{n}(:,3) = Qdls;
end
%% Plot out trial max
dat1 = QM{3}(:,3);
dat1 = rmoutliers(dat1);
dat2 = rmoutliers(dat2);
dat2 = QD{5}(:,3);
temp =nan(max([length(dat1),length(dat2)]),2);
temp(1:length(dat1),1)= dat1;
temp(1:length(dat2),2) = dat2;
temp = rmoutliers(temp);
figure,customBarplot(temp);
box off,axis square, set(gca,'tickdir','out')
xlim([0.0 3])
%% Plot it out
cellCount = numel(QM);                % Number of cells
pulls = size(QM{1}, 2);               % Number of pulls

pullMeans = nan(cellCount, pulls);              % Preallocate for cell means

for i = 1:cellCount
    pullMeans(i, :) = mean(QM{i}, 1); % Mean across trials
end

mean_vals = mean(pullMeans, 1);
sem_vals = std(pullMeans, 0, 1) / sqrt(cellCount);

figure;
hold on;
colors = lines(cellCount);
% Plot individual cell lines and dots
for i = 1:cellCount
    plot(1:pulls, pullMeans(i, :), '-o', ...
        'Color', colors(i, :), ... 
        'MarkerFaceColor', colors(i, :), ...
        'MarkerEdgeColor', colors(i, :), ...
        'LineWidth', 1.5, 'MarkerSize', 7);
end
% Overlay the barplot (mean ± SEM)
bar_handle = bar(1:pulls, mean_vals, 0.7, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
for b = 1:pulls
    bar_handle.FaceColor = [0.4 0.4 0.4]; % Optional: neutral bar color
end
xlabel('Pull');
ylabel('Average Modulation');
set(gca, 'XTick', 1:pulls,'tickdir','out');
hold off;
xlim([0 4])
axis square
%% Plot it out
cellCount = numel(QD);                % Number of cells
pulls = size(QD{1}, 2);               % Number of pulls

pullMeans = nan(cellCount, pulls);              % Preallocate for cell means

for i = 1:cellCount
    pullMeans(i, :) = mean(QD{i}, 1); % Mean across trials
end

mean_vals = mean(pullMeans, 1);
sem_vals = std(pullMeans, 0, 1) / sqrt(cellCount);

figure;
hold on;
colors = lines(cellCount);
% Plot individual cell lines and dots
for i = 1:cellCount
    plot(1:pulls, pullMeans(i, :), '-o', ...
        'Color', colors(i, :), ... 
        'MarkerFaceColor', colors(i, :), ...
        'MarkerEdgeColor', colors(i, :), ...
        'LineWidth', 1.5, 'MarkerSize', 7);
end
% Overlay the barplot (mean ± SEM)
bar_handle = bar(1:pulls, mean_vals, 0.7, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
for b = 1:pulls
    bar_handle.FaceColor = [0.4 0.4 0.4]; % Optional: neutral bar color
end
xlabel('Pull');
ylabel('Average Modulation');
set(gca, 'XTick', 1:pulls,'tickdir','out');
hold off;
xlim([0 4])
axis square
%% LOCAL function
function modulation = getModulation(neuralDynamics)
modulation = [];
pullSpeed = squeeze(neuralDynamics.warpedSpikes.GPFA.p1neuralDynamics.hitOnly.speed.speed(1,:,:));
time = neuralDynamics.warpedSpikes.warpedTime(:,1);
modulation(:,1) = getTrajModulation(pullSpeed,time);

pullSpeed = squeeze(neuralDynamics.warpedSpikes.GPFA.p2neuralDynamics.hitOnly.speed.speed(1,:,:));
time = neuralDynamics.warpedSpikes.warpedTime(:,2);
modulation(:,2) = getTrajModulation(pullSpeed,time);

pullSpeed = squeeze(neuralDynamics.warpedSpikes.GPFA.p3neuralDynamics.hitOnly.speed.speed(1,:,:));
time = neuralDynamics.warpedSpikes.warpedTime(:,3);
modulation(:,3) = getTrajModulation(pullSpeed,time);

figure,
customBarplot(modulation)
ylabel('Speed Modulation'),axis square
set(gca,'tickdir','out'),box off
end

function modulation = getTrajModulation(speed,time)
time = linspace(time(1),time(end),size(speed,1));
bWin = time<-1.5; % This is the fixed baseline because we want the initial condition prior to pull. 
baseline = mean(speed(bWin));
speed = mean(speed(time>-1.5 & time<1.5,:))-baseline;
% [~,idx] = max(abs(speed),[],1);
% n = size(speed,2);                   % number of columns (trials)
% vals = speed(sub2ind(size(speed), idx, 1:n));  % extract true values (with sign)
modulation = speed/baseline;  
end

function prepRSLDS(warpedSpikes,filename)
[fpath,fname,ext] = fileparts(filename)
dat = warpedSpikes;
filename = fullfile(fpath,[fname,'_rslds.npy']);
writeNPY(dat, filename)
disp('rslds data saved')
end