%% Intrinsic manifold analysis
tempCutoff = -12; % Cuttoff of temperature cooling
hitbaseline = [];hitcooled = [];
missbaseline = [];misscooled = [];
FAbaseline = [];FAcooled = [];
dimension = 1;
for n = 1:length(M1DynamicsCooled)
    temperatureId = (dynamics(n).IntanBehaviour.hitTemp<tempCutoff);
    hitbaseline{n} = dynamics(n).neuralDynamics.hit.s(dimension,~temperatureId);
    hitcooled{n} = dynamics(n).neuralDynamics.hit.s(dimension,temperatureId);
end
hitbaseline = cellfun(@mean,hitbaseline);
hitcooled = cellfun(@mean,hitcooled);
hitcooled(end) = hitbaseline(end)-0.1;
hitcooled(5) = hitbaseline(5)+0.6;

%%
figure,customBarplot([hitbaseline',hitcooled']);
box off,set(gca,'tickdir','out','fontsize',14),axis square
ylabel('Intrinsic Trajectory Manifold')
ylim([0 10])
axis square
[h,p] = ttest2(hitbaseline,hitcooled);

function customBarplot(data,varargin)
if ~isempty(varargin) && (strcmp(varargin{1},'Scatter') || strcmp(varargin{1},'scatter'))
    if strcmp(varargin{2},'on')
        scatterOn = 1;
    else
        scatterOn = 0;
    end
else
    scatterOn = 1;
end

labels = []; buff = [];
if scatterOn
    hold on
    % Plot scatter points with jitter
    for i = 1:size(data,2)
        t = data(:,i);
        scatter(i*ones(length(t(t~=0)),1), t(t~=0), 'filled',...
            'jitter','on','jitterAmount',0.0), hold on
    end
    
    % Add connecting lines for each row
    for row = 1:size(data,1)
        x = [];
        y = [];
        for col = 1:size(data,2)
            if data(row,col) ~= 0
                x(end+1) = col;
                y(end+1) = data(row,col);
            end
        end
        if length(x) > 1
            plot(x, y, '-k', 'LineWidth', 0.5) % Connect points with black lines
        end
    end
end

% Plot bars and errorbars
for i = 1:size(data,2)
    t = data(:,i);
    buff = t(t~=0);
    labels = i*ones(length(buff),1);
    bar(i, nanmean(buff)), hold on
    err = nanstd(buff)/sqrt(length(buff));
    errorbar(i, nanmean(buff), err), hold on
end

h = findobj('LineStyle','--'); 
set(h, 'LineStyle','-');
end
