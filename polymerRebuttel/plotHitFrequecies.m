%% Main Analysis Script
% Extract and preprocess data
leverTracePre = arrayfun(@(x) x.rawtrace, BehaviourPreLaser.hitTrace,'UniformOutput',false);
leverTracePre = abs(horzcat(leverTracePre{:})-133);

leverTraceLaser = arrayfun(@(x) x.rawtrace, BehaviourLaser.hitTrace,'UniformOutput',false);
leverTraceLaser = abs(horzcat(leverTraceLaser{:})-133);

leverTracePost = arrayfun(@(x) x.rawtrace, BehaviourPostLaser.hitTrace,'UniformOutput',false);
leverTracePost = abs(horzcat(leverTracePost{:})-133);

% Smooth traces
smoothedPre   = smoothdata(leverTracePre, 2, 'gaussian', 10);
smoothedLaser = smoothdata(leverTraceLaser, 2, 'gaussian', 10);
smoothedPost  = smoothdata(leverTracePost, 2, 'gaussian', 10);

% Calculate metrics (Peak Speed and AUC)
peakSpeedPre = max(diff(smoothedPre,1,2), [], 2);
aucPre = trapz(smoothedPre, 2);

peakSpeedLaser = max(diff(smoothedLaser,1,2), [], 2);
aucLaser = trapz(smoothedLaser, 2);

peakSpeedPost = max(diff(smoothedPost,1,2), [], 2);
aucPost = trapz(smoothedPost, 2);

%% Create Integrated Figure
figure('Color','w','Position',[100 100 1200 800]);
tiledlayout(2,3,'TileSpacing','compact');

% Individual Traces
nexttile(1)
plot(smoothedPre', 'Color',[0.5 0.5 0.5 0.25])
title('Pre-Laser Traces'), axis tight

nexttile(2)
plot(smoothedLaser', 'Color',[0.9 0.4 0.2 0.25])
title('LaserOn Traces'), axis tight

nexttile(3)
plot(smoothedPost', 'Color',[0.2 0.8 0.3 0.25])
title('Post-Laser Traces'), axis tight

% Mean Traces
nexttile(4,[1 3])
hold on
plot(mean(smoothedPre,1), 'Color',[0.3 0.3 0.3], 'LineWidth',2)
plot(mean(smoothedLaser,1), 'Color',[0.9 0.4 0.2], 'LineWidth',2)
plot(mean(smoothedPost,1), 'Color',[0.2 0.8 0.3], 'LineWidth',2)
legend({'PreLaser','LaserOn','PostLaser'},'Location','northeastoutside')
title('Mean Traces'), grid on, axis tight

% Violin Plots with Statistics
for i = 5:6
    nexttile(i)
    if i == 5
        data = {peakSpeedPre, peakSpeedLaser, peakSpeedPost};
        title('Peak Pull Speed')
    else
        data = {aucPre, aucLaser, aucPost};
        title('Area Under Curve (AUC)')
    end
    
    % Combine data
    allData = [data{1}; data{2}; data{3}];
    groups = [ones(numel(data{1}),1); 2*ones(numel(data{2}),1); 3*ones(numel(data{3}),1)];
    
    % Plot violins
    violinplot(allData, categorical(groups,[1 2 3],{'PreLaser','LaserOn','PostLaser'}), ...
        'ViolinColor',[0.2 0.6 0.8; 0.9 0.4 0.2; 0.6 0.8 0.3]);
    
    % Add statistics
    [~,~,stats] = kruskalwallis(allData, groups, 'off');
    [results,~,~] = multcompare(stats, 'Display','off');
    y_max = max(allData)*1.1;
    
    % Add significance bars
    hold on
    for r = 1:size(results,1)
        if results(r,6) < 0.05
            plot(results(r,1:2), [1 1]*y_max*(1+r*0.05), 'k-')
            text(mean(results(r,1:2)), y_max*(1.05+r*0.05), ...
                getSigSymbol(results(r,6)), 'HorizontalAlignment','center')
        end
    end
    
    % Add sample sizes
    text(1, y_max, ['N=' num2str(numel(data{1}))], 'HorizontalAlignment','center')
    text(2, y_max, ['N=' num2str(numel(data{2}))], 'HorizontalAlignment','center')
    text(3, y_max, ['N=' num2str(numel(data{3}))], 'HorizontalAlignment','center')
    axis square
end

%% Helper for significance symbols (inline)
function symbol = getSigSymbol(p)
    if p < 0.001
        symbol = '***';
    elseif p < 0.01
        symbol = '**';
    elseif p < 0.05
        symbol = '*';
    else
        symbol = '';
    end
end
