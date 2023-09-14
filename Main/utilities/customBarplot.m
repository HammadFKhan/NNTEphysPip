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

labels = [];buff = [];
for i = 1:size(data,2)
    t = data(:,i);
    buff = t(t~=0);
    labels = i*ones(length(buff),1);
    bar(i,mean(buff)), hold on
    err = std(buff)/sqrt(length(buff));
    errorbar(i,mean(buff),err),hold on
end

h=findobj('LineStyle','--'); set(h, 'LineStyle','-');

if scatterOn
    for i = 1:size(data,2)
        t = data(:,i);
        scatter(i*ones(length(t(t~=0)),1),t(t~=0),'filled','jitter','on','jitterAmount',0.1)
    end
end
