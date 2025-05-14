%% Polymer rebuttel

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
fname = 'Y:\Om\Behaviour\PolymerRebuttel\polymerMouse22025_05_11_04.27.PM.csv'; %% Control
fname = 'Y:\Om\Behaviour\PolymerRebuttel\polymerMouse22025_05_12_02.07.PM.csv'; %% Polymer
[Behaviour] = readLever(parameters,[],fname);
%% Seperate trials into laser on and laser off trials
BehaviourPreLaser = Behaviour;
BehaviourLaser = Behaviour;
BehaviourPostLaser = Behaviour;

laseridx = find(abs(diff(Behaviour.B(:,5)))==1);

BehaviourPreLaser.hitTrace(arrayfun(@(x) x.i1>laseridx(1), BehaviourPreLaser.hitTrace)) = []; % Remove trials afer laser is turned on
BehaviourLaser.hitTrace(arrayfun(@(x) x.i1<laseridx(1), BehaviourLaser.hitTrace)) = []; % Remove trials before laser is turned on and after it is off
BehaviourLaser.hitTrace(arrayfun(@(x) x.i1>laseridx(2), BehaviourLaser.hitTrace)) = []; % Remove trials before laser is turned on and after it is off
BehaviourPostLaser.hitTrace(arrayfun(@(x) x.i1<laseridx(2), BehaviourPostLaser.hitTrace)) = []; % Remove trials before laser is turned off

BehaviourPreLaser.missTrace(arrayfun(@(x) x.i1>laseridx(1), BehaviourPreLaser.missTrace)) = []; % Remove trials afer laser is turned on
BehaviourLaser.missTrace(arrayfun(@(x) x.i1<laseridx(1), BehaviourLaser.missTrace)) = []; % Remove trials before laser is turned on and after it is off
BehaviourLaser.missTrace(arrayfun(@(x) x.i1>laseridx(2), BehaviourLaser.missTrace)) = []; % Remove trials before laser is turned on and after it is off
BehaviourPostLaser.missTrace(arrayfun(@(x) x.i1<laseridx(2), BehaviourPostLaser.missTrace)) = []; % Remove trials before laser is turned off
%% Hit rate and lever response
hitTimePre = arrayfun(@(x) x.t1, BehaviourPreLaser.hitTrace);
hitTimePre = diff(hitTimePre);
hitTimePre(abs(hitTimePre)>100) = [];

hitTimeLaser = arrayfun(@(x) x.t1, BehaviourLaser.hitTrace);
hitTimeLaser = diff(hitTimeLaser);
hitTimeLaser(abs(hitTimeLaser)>100) = [];

hitTimePost = arrayfun(@(x) x.t1, BehaviourPostLaser.hitTrace);
hitTimePost = diff(hitTimePost);
hitTimePost(abs(hitTimePost)>100) = [];

missTimePre = arrayfun(@(x) x.t1, BehaviourPreLaser.missTrace);
missTimePre = diff(missTimePre);
missTimePre(abs(missTimePre)>100) = [];

missTimeLaser = arrayfun(@(x) x.t1, BehaviourLaser.missTrace);
missTimeLaser = diff(missTimeLaser);
missTimeLaser(abs(missTimeLaser)>100) = [];

missTimePost = arrayfun(@(x) x.t1, BehaviourPostLaser.missTrace);
missTimePost = diff(missTimePost);
missTimePost(abs(missTimePost)>100) = [];

getTaskRate(hitTimePre,hitTimeLaser,hitTimePost)

getTaskRate(missTimePre,missTimeLaser,missTimePost)

%% LOCAL FUNCTION
function getTaskRate(hitTimePre,hitTimeLaser,hitTimePost)
% Your existing code to extract and clean hit times

% Organize data for statistical testing
allHitTimes = [hitTimePre, hitTimeLaser, hitTimePost];
groupLabels = [ones(length(hitTimePre), 1); 
               2*ones(length(hitTimeLaser), 1);
               3*ones(length(hitTimePost), 1)];

% Perform Kruskal-Wallis test (non-parametric alternative to one-way ANOVA)
[p_kw, tbl_kw, stats_kw] = kruskalwallis(allHitTimes, groupLabels, 'off');

% Perform multiple comparisons to determine which groups differ
[results, ~, ~, gnames] = multcompare(stats_kw, 'Display', 'off');

% Create a beautiful violin plot
figure('Position', [100, 100, 800, 600], 'Color', 'white');

% Create categorical array for labels
groupLabels = categorical(groupLabels, [1,2,3], {'PreLaser', 'LaserOn', 'PostLaser'});

% Plot with enhanced aesthetics
violinplot(allHitTimes, groupLabels, ...
    'ShowData', true, ...                        
    'ViolinColor', [0.2 0.6 0.8; 0.9 0.4 0.2; 0.6 0.8 0.3], ... 
    'ViolinAlpha', 0.25, ...                      
    'ShowBox', false, ...                        
    'ShowMean', true, ...                        
    'ShowMedian', true);                        

% Improve title and labels
title('Laser Blanking Control', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Perturbation Phase', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Intertrial Hit Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% Add significance bars
hold on;
y_max = max(allHitTimes) * 1.1;
y_positions = [y_max*1.05, y_max*1.1, y_max*1.15];
group_positions = [1, 2, 3]; % x-positions of groups

% Add significance bars and p-values between groups
for i = 1:size(results,1)
    group1 = results(i,1);
    group2 = results(i,2);
    p_val = results(i,6);
    
    % Draw significance bars if significant
    if p_val < 0.05
        x1 = group_positions(group1);
        x2 = group_positions(group2);
        y = y_positions(i);
        
        % Plot the bar
        plot([x1, x2], [y, y], 'k-', 'LineWidth', 1.5);
        plot([x1, x1], [y-y_max*0.01, y], 'k-', 'LineWidth', 1.5);
        plot([x2, x2], [y-y_max*0.01, y], 'k-', 'LineWidth', 1.5);
        
        % Add significance markers
        if p_val < 0.001
            text((x1+x2)/2, y+y_max*0.01, '***', 'HorizontalAlignment', 'center', 'FontSize', 12);
        elseif p_val < 0.01
            text((x1+x2)/2, y+y_max*0.01, '**', 'HorizontalAlignment', 'center', 'FontSize', 12);
        else
            text((x1+x2)/2, y+y_max*0.01, '*', 'HorizontalAlignment', 'center', 'FontSize', 12);
        end
    end
end

% Add Kruskal-Wallis test result
text(1, y_max*1.2, ['Kruskal-Wallis test: p = ' num2str(p_kw, '%.4f')], 'FontSize', 10);

% Adjust y-axis to accommodate significance bars
curr_ylim = ylim;
ylim([curr_ylim(1), y_positions(end)+y_max*0.1]);

% Add box around the plot
box on;
set(gca, 'LineWidth', 1.5);
y_max = max(allHitTimes) * 1.1;
Ns = [length(hitTimePre), length(hitTimeLaser), length(hitTimePost)];
for i = 1:3
    text(i, y_max, ['N = ' num2str(Ns(i))], ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 12, 'FontWeight', 'bold');
end

axis square
hold off;
end
