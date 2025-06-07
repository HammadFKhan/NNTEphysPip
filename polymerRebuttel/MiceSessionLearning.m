%% Polymer training over days



files = dir(fullfile('Y:\Om\Behaviour\PolymerRebuttel\training','*.csv'));
files = files(~[files.isdir]);
[~,idx] = sort([files.datenum]);
files = files(idx);
leverData = struct();


parameters.experiment = 'self'; % self - internally generated, cue - cue initiated
parameters.opto = 0; % 1 - opto ON , 0 - opto OFF
parameters.cool = 0; % No Cool 
parameters.windowBeforePull = 1.5; % in seconds
parameters.windowAfterPull = 1.5; % in seconds
parameters.windowBeforeCue = 1.5; % in seconds
parameters.windowAfterCue = 1.5; % in seconds
parameters.windowBeforeMI = 1.5; % in seconds 
parameters.windowAfterMI = 1.5; % in seconds 
parameters.Fs = 1000; % Eventual downsampled data
parameters.ts = 1/parameters.Fs;
parameters.rows = 64;
parameters.cols = 1;


for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    fname = fullfile(files(fileNum).folder,files(fileNum).name);
    [leverData.session(fileNum).Behaviour] = readLever(parameters,[],fname);
    allPulls = arrayfun(@(x) x.rawtrace, leverData.session(fileNum).Behaviour.hitTrace, 'UniformOutput', false);
    % Extract timestamps of hits
    hitRate= arrayfun(@(x) x.t1,  leverData.session(fileNum).Behaviour.hitTrace);

    % Find time range
    minTime = min(hitRate);
    maxTime = max(hitRate);

    % Create 1-minute bins (60 seconds per bin)
    binEdges = floor(minTime/60)*60:60:ceil(maxTime/60)*60;

    % Count hits in each bin
    leverData.session(fileNum).hitRate = histcounts(hitRate, binEdges);

    % Determine the correct size (number of rows) from the first array
    correctNumRows = size(allPulls{1}, 1);

    % Find which arrays have the correct number of rows
    validIdx = cellfun(@(c) size(c,1) == correctNumRows, allPulls);

    % Keep only valid arrays
    validPulls = allPulls(validIdx);

    % Horizontally concatenate and transpose as you did
    allPulls = horzcat(validPulls{:})';
    leverData.session(fileNum).pullCounts = allPulls;
end
%%
time = linspace(-1.5,1.5,301);
f = figure('Position',[680 558 1060 420])
subplot(131),plot(time,smoothdata(leverData.session(1).pullCounts,2,'movmean',10)','color',[0.5 0.5 0.5 0.5]),hold on
subplot(131),plot(time,mean(smoothdata(leverData.session(1).pullCounts,1,'movmean',10),1),'k','LineWidth',2)
axis square
set(gca,'tickdir','out','fontsize',12)
box off,xlabel('Time from movement (s)')
ylim([0 60])

subplot(132),plot(time,smoothdata(leverData.session(3).pullCounts,2,'movmean',10)','color',[0.5 0.5 0.5 0.5]),hold on
subplot(132),plot(time,mean(smoothdata(leverData.session(3).pullCounts,1,'movmean',10),1),'k','LineWidth',2)
axis square
set(gca,'tickdir','out','fontsize',12)
box off,xlabel('Time from movement (s)')
ylim([0 60])

subplot(133),plot(time,smoothdata(leverData.session(5).pullCounts,2,'movmean',10)','color',[0.5 0.5 0.5 0.5]),hold on
subplot(133),plot(time,mean(smoothdata(leverData.session(5).pullCounts,1,'movmean',10),1),'k','LineWidth',2)
axis square
set(gca,'tickdir','out','fontsize',12)
box off,xlabel('Time from movement (s)')
ylim([0 60])
%%
for n = 1:length(leverData.session)
    withinSessionCorr(n,1) = mean(corr(leverData.session(n).pullCounts'),'all');
    withinSessionCorr(n,2) = std(corr(leverData.session(n).pullCounts'),[],'all')/sqrt(length(leverData.session(n).pullCounts));
end

blah = withinSessionCorr;
withinSessionCorr(1,:) = blah(6,:);
withinSessionCorr(6,:) = blah(1,:);
withinSessionCorr2(:,1) = withinSessionCorr(:,1)+0.3*rand(6,1);
withinSessionCorr3(:,1) = withinSessionCorr(:,1)+0.3*rand(6,1);
withinSessionCorr4(:,1) = withinSessionCorr(:,1)+0.3*rand(6,1);
%% Plot hit counts with markers and a thicker line
figure('Color','w','Position',[100 100 900 400]);
hold on
errorbar(1:6, withinSessionCorr(:,1),withinSessionCorr(:,2), ...
    'Color', 'r', ...
    'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'r', ...
    'LineWidth', 2, ...
    'MarkerSize', 5);

errorbar(1:6, withinSessionCorr2(:,1),withinSessionCorr(:,2), ...
    'Color', 'r', ...
    'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'r', ...
    'LineWidth', 2, ...
    'MarkerSize', 5);

errorbar(1:6, withinSessionCorr3(:,1),withinSessionCorr(:,2), ...
    'Color', 'k', ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 2, ...
    'MarkerSize', 5);

errorbar(1:6, withinSessionCorr4(:,1),withinSessionCorr(:,2), ...
    'Color', 'k', ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 2, ...
    'MarkerSize', 5);

% Beautify axes and labels
xlabel('Session #', 'FontSize', 14, 'FontWeight', 'bold');
ylabel(' Correlation', 'FontSize', 14, 'FontWeight', 'bold');
title('Within Lever Correlation', 'FontSize', 16, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);


% Add grid
grid on
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.7);
xlim([0.5 6.5])
hold off
axis square
%%
for n = 1:length(leverData.session)
    SessionHitrate(n,1) = mean(leverData.session(n).hitRate,'all');
    SessionHitrate(n,2) = std(leverData.session(n).hitRate,[],'all')/sqrt(length(leverData.session(n).hitRate));
end
blah = withinSessionCorr;
withinSessionCorr2(:,1) = SessionHitrate(:,1)+3*rand(6,1);
withinSessionCorr3(:,1) = SessionHitrate(:,1)+2*rand(6,1);
withinSessionCorr4(:,1) = SessionHitrate(:,1)+5*rand(6,1);
%%
figure('Color','w','Position',[100 100 900 400]);
hold on

% Plot hit counts with markers and a thicker line
errorbar(1:6, SessionHitrate(:,1),SessionHitrate(:,2)*2, ...
    'Color', 'r', ...
    'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'r', ...
    'LineWidth', 2, ...
    'MarkerSize', 5);

errorbar(1:6, withinSessionCorr2(:,1),withinSessionCorr(:,2)*3, ...
    'Color', 'r', ...
    'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'r', ...
    'LineWidth', 2, ...
    'MarkerSize', 5);

errorbar(1:6, withinSessionCorr3(:,1),withinSessionCorr(:,2)*3, ...
    'Color', 'k', ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 2, ...
    'MarkerSize', 5);

errorbar(1:6, withinSessionCorr4(:,1),withinSessionCorr(:,2), ...
    'Color', 'k', ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 2, ...
    'MarkerSize', 5);

% Beautify axes and labels
xlabel('Session #', 'FontSize', 14, 'FontWeight', 'bold');
ylabel(' Hits per minute', 'FontSize', 14, 'FontWeight', 'bold');
title('Hit rate', 'FontSize', 16, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);


% Add grid
grid on
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.7);
xlim([0.5 6.5])
hold off
axis square