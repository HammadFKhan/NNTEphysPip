function spikeRaster(Spikes)
%Create Spike rasters of units by depth during brhavior state
spikeTemp = squeeze(struct2cell(Spikes.spikes.L23Run)); % Spike temp for layer 23 of cell cell array
spikeTemp = vertcat(spikeTemp{:});
figure,hold on
for i = 1:size(spikeTemp,1)
    r = horzcat(spikeTemp(i,:));
    raster{i} = vertcat(r{:});
    plot(raster{i},(i/5)*ones(size(raster{i},1)),'.')
end

spikeTemp = squeeze(struct2cell(Spikes.spikes.L4Run)); % Spike temp for layer 23 of cell cell array
spikeTemp = vertcat(spikeTemp{:});
figure,hold on
for i = 1:size(spikeTemp,1)
    r = horzcat(spikeTemp(i,:));
    raster{i} = vertcat(r{:});
    plot(raster{i},(i/5)*ones(size(raster{i},1)),'.')
end
spikeTemp = [];

spikeTemp = squeeze(struct2cell(Spikes.spikes.L5Run)); % Spike temp for layer 23 of cell cell array
spikeTemp = vertcat(spikeTemp{:});
figure,hold on
for i = 1:size(spikeTemp,1)
    r = horzcat(spikeTemp(i,:));
    raster{i} = vertcat(r{:});
    plot(raster{i},(i/5)*ones(size(raster{i},1)),'.')
end
spikeTemp = [];
end
