function spikeImage = spike_map(DeltaFoverF,time,color)
if nargin < 3 || strcmp(color,'')
    color = jet;
end
xyv = DeltaFoverF;
x = length(xyv(1,:)); % Extract x value
y = length(xyv(:, 1)); % Extract y value
v = xyv; % Extract intensity value
minX = min(x);
maxX = max(x);
minY = min(y);
maxY = max(y);
minV = min(v,[],'all');
maxV = max(v,[],'all');
avgV = mean(v');
spikeImage = zeros(y,x);
for i = 1:x
    for ii = 1:y
        spikeImage(ii,i) = DeltaFoverF(ii,i);
    end
end

% if y<25
ax = time;
[grad,~]=colorGradient([.1 .1 .1],[1 1 1],64);
colormap(jet);
imagesc(spikeImage,'Xdata',time);
% h = colorbar;
% set(get(h,'title'),'string','FR');
axis on;axis tight;box off;
xlabel('Time (s)'); 
ylabel('Neuron');
% else
%     imshow(spikeImage);axis on;axis tight;box off;
%     h = colorbar;
%     set(get(h,'title'),'string','\Delta F/F (%)');
%     xlabel('Time (s)');
%     ylabel('Neuron');
% end

% Colormap is not gray scale.
% Apply some other colormap if you want
end
