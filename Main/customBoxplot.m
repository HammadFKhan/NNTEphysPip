function customBoxplot(data)
%% Custom boxplot that plots with overlayed scatter plot

% Parse input
if length(size(data))>1
    labels = [];buff = [];
    for i = 1:size(data,2) 
        t = data(:,i);
        buff = [buff;t(t>0)];
        labels = [labels;repmat({num2str(i)},length(t(t>0)),1)];
    end
    boxplot(buff,labels), hold on
    for i = 1:size(data,2)
        t = data(:,i);
        buff = [buff;t(t>0)];
        labels = repmat({num2str(i)},length(buff),1);
        scatter(i*ones(length(t(t>0)),1),t(t>0),'filled','jitter','on','jitterAmount',0.08)
    end
else
    boxplot(data),hold on
    scatter(ones(size(data,1),1),data,'filled','jitter','on','jitterAmount',0.08)
end
