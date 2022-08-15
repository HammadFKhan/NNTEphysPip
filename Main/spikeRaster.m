function spikeRaster(Spikes)
%Create Spike rasters of units by depth during brhavior state
spikeTemp = squeeze(struct2cell(Spikes.spikes.L23Rest)); % Spike temp for layer 23 of cell cell array
spikeTemp = vertcat(spikeTemp{:});
figure,hold on
for i = 1:11
    r = horzcat(spikeTemp(i,:));
    idx = cellfun('size',r,1);
    r(cellfun('size',r,1)>10) = [];
    raster = vertcat(r{:});
    scatter(raster,i*ones(size(raster,1),1),'filled','r')
    rasterAvg{i} = raster;
end

% spikeTemp = squeeze(struct2cell(Spikes.spikes.L4Run)); % Spike temp for layer 23 of cell cell array
% spikeTemp = vertcat(spikeTemp{:});
% figure,hold on
% for i = 1:size(spikeTemp,1)
%     r = horzcat(spikeTemp(i,:));
%     raster{i} = vertcat(r{:});
%     plot(raster{i},(i/5)*ones(size(raster{i},1)),'.')
% end
% spikeTemp = [];

spikeTemp = squeeze(struct2cell(Spikes.spikes.L5Rest)); % Spike temp for layer 23 of cell cell array
spikeTemp = vertcat(spikeTemp{:});
figure,hold on
for i = 1:11
    r = horzcat(spikeTemp(i,:));
    idx = cellfun('size',r,1);
    r(cellfun('size',r,1)>20) = [];
    raster = vertcat(r{:});
    scatter(raster,i*ones(size(raster,1),1),'filled','b')
end
spikeTemp = [];
end
