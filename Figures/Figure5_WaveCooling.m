M1WavesCoolingFiles = readcell("\\10.165.57.13\Sutter_backup\Hammad\Ephys\LeverTask\Data_for_Figures\WavesCooling\WavesCoolingFiles.xlsx");
M1WavesCoolingFiles(1,:) = []; % Delete the first row 
M1RTCooling = [];
% commenttorunsafetyswitch
%% Reading files from each locations
fileNo = 1;
for i=1:size(M1WavesCoolingFiles,1)
    if(M1WavesCoolingFiles{i,3} == 0)
        continue;
    end
    fileName = [M1WavesCoolingFiles{i,2} '\' M1WavesCoolingFiles{i,1}];
    M1RTCooling(fileNo).fileID = i;
    M1RTCooling(fileNo).fileName = fileName;

    disp(['Reading from ' fileName]);
    matObj = matfile(fileName);
    variables = who(matObj);
    skip = 0;
    % Getting IntanBehaviour
    if ismember('IntanBehaviour', variables)
        disp("Loading variables - IntanBehaviour and parameters");
        parameters = load(fileName,"parameters");
        parameters = parameters.parameters;
        IntanBehaviour = load(fileName,"IntanBehaviour");
        IntanBehaviour = IntanBehaviour.IntanBehaviour;
    else
        disp("Could not find IntanBehaviour");
        IntanBehaviour = [];
        skip = 1;
    end
    M1RTCooling(fileNo).parameters = parameters;
    baselineTemp = 20;coolingTemp = 14;recoveryTemp = 18;
    [IntanBehaviourBaseline, ~, ~, ~, ~, ~] = splitCooling(baselineTemp, coolingTemp, recoveryTemp, IntanBehaviour,[]);
    
    bodyTempTime = 2*60; % in seconds
    bodyTemp = mean(IntanBehaviour.tempTrace(1:bodyTempTime*parameters.Fs),'all');
    hitTemp = horzcat(IntanBehaviour.cueHitTrace.temp);
    RT = horzcat(IntanBehaviour.cueHitTrace.reactionTime);
    muRTBaseline = mean(horzcat(IntanBehaviourBaseline.cueHitTrace.reactionTime),'all');
    stdRTBaseline = std(horzcat(IntanBehaviourBaseline.cueHitTrace.reactionTime),0,'all');
    M1RTCooling(fileNo).zRT = (RT-muRTBaseline)/stdRTBaseline;   
    M1RTCooling(fileNo).delTemp = hitTemp-bodyTemp;
    M1RTCooling(fileNo).Temp = hitTemp;
    M1RTCooling(fileNo).RT = RT;
    fileNo = fileNo + 1;
    clear IntanBehaviour IntanBehaviourBaseline;
end


%% Combining all RT and delTemp
load("\\10.165.57.13\Sutter_backup\Hammad\Ephys\LeverTask\Data_for_Figures\WavesCooling\M1RTCoolingCombinedNew.mat")
RT = horzcat(M1RTCooling.RT);
zRT = horzcat(M1RTCooling.zRT);
delTemp = horzcat(M1RTCooling.delTemp);

% Rejecting any glitches in RT 
k = find(RT<0);
zRT(k) = [];delTemp(k) = [];

% Rejecting outliers in RT (99.99 cutoff)
k = find(zRT<-3.891 | zRT>3.891);
zRT(k) = [];delTemp(k) = [];

mdl = fitlm(delTemp,zRT)

% Define custom blue and red as RGB
rgb_red = [190 31 89]/255;    % #be1f59
rgb_blue = [40 153 196]/255;  % #2899c4

% Number of colors for the gradient
nColor = 256;

% Interpolate colormap from blue to red
cmap = [linspace(rgb_blue(1), rgb_red(1), nColor)', ...
        linspace(rgb_blue(2), rgb_red(2), nColor)', ...
        linspace(rgb_blue(3), rgb_red(3), nColor)'];

% Normalize temperature data
delTempNorm = (delTemp - min(delTemp)) / (max(delTemp) - min(delTemp));
colorIdx = max(1, min(nColor, round(1 + delTempNorm * (nColor-1)))); % Ensure in [1 nColor]

figure; hold on;

% Scatter plot with colored dots
scatter(delTemp, zRT, 40, cmap(colorIdx,:), 'filled');

% Plot regression line
xfit = linspace(min(delTemp), max(delTemp), 100);
yfit = predict(mdl, xfit');
plot(xfit, yfit, 'k-', 'LineWidth', 2);

xlabel('Temperature Change');
ylabel('Reaction Time (z-scored)');
box off; set(gca, 'TickDir', 'out', 'FontSize', 14);

% Add custom colorbar
colormap(cmap);
cb = colorbar;
cb.Label.String = 'Temperature Change';
caxis([min(delTemp) max(delTemp)]);
xlim([-17 1]); ylim([-3 4])
hold off;box off;
%%
figure,histogram(delTemp,-17:1,'Normalization','probability','LineStyle','none'),axis square
xlim([-17 1]);box off,set(gca,'tickdir','out')
