for i = 1:1000
    x = 0:0.05:10;
    total(i,:) = 3*sin(x+50*rand(1));
    total2(i,:) = 1*sin(x+5*rand(1));
end
total = [total;total2];
%%
total = [];
chunk = 1:floor(1024/2):size(linearProbe,2);
train_labels = {};
total = zeros(floor(1024/2),length(chunk)-1,2);
count = 1;
for electrode = [4 20]
    disp(['Analyzing Electrode: ' num2str(electrode)]);
    for i = 1:length(chunk)-1
        total(1:floor(1024/2),i,count) = linearProbe(electrode,chunk(i):chunk(i+1)-1);
    end
    train_labels = [train_labels;repmat({num2str(electrode)},length(chunk)-1,1)];
    count = count+1;
end

total = reshape(total,floor(1024/2),[]);
%%
total2 = [];
for i = 1:size(total,2);
    wavelet = cwt(total(:,1),1024,'FrequencyLimit',[5 80]);
    total2 = [total2;abs(wavelet)];
end
%%
no_dims = 2;
initial_dims = 10;
perplexity = 10;

% Run t-SNE
mappedX = tsne(total,...
    'Algorithm','barneshut',...
    'NumDimensions',no_dims,...
    'NumPCAComponents',initial_dims,...
    'Perplexity', perplexity,...
    'Exaggeration',4,...
    'Verbose',2);
figure,gscatter(mappedX(:,1), mappedX(:,2));
%%
colors = hsv(22);
figure,
for i = 1:22
    idx = find(strcmp(train_labels,num2str(i)));
    subplot(6,4,i),scatter(mappedX(idx,1),mappedX(idx,2),3,'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',colors(i,:))
    axis([-80 80 -100 100])
end

%%
train_labels = cell(1024,1);
for i = 1:1024
    train_labels{i} = num2str(i);
end