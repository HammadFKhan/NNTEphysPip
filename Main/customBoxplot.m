function customBoxplot(data,varargin)
%% Custom boxplot that plots with overlayed scatter plot
if ~isempty(varargin) && (strcmp(varargin{1},'Scatter') || strcmp(varargin{1},'scatter'))
    if strcmp(varargin{2},'on')
        scatterOn = 1;
    else
        scatterOn = 0;
    end
else
    scatterOn = 1;
end
% Parse input
if length(size(data))>1
    labels = [];buff = [];
    for i = 1:size(data,2)
        t = data(:,i);
        buff = [buff;t(t>0)];
        labels = [labels;repmat({num2str(i)},length(t(t>0)),1)];
    end
    boxplot(buff,labels), hold on
    h=findobj('LineStyle','--'); set(h, 'LineStyle','-');

    if scatterOn
        for i = 1:size(data,2)
            t = data(:,i);
            buff = [buff;t(t>0)];
            labels = repmat({num2str(i)},length(buff),1);
            scatter(i*ones(length(t(t>0)),1),t(t>0),'filled','jitter','on','jitterAmount',0.1)
        end
    end
else
    boxplot(data),hold on
    if scatterOn
        scatter(ones(size(data,1),1),data,'filled','jitter','on','jitterAmount',0.08)
    end
end

