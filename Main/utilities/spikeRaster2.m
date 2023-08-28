function [L23rasterAvg,L5rasterAvg] = spikeRaster2(Spikes, flag,sz)
%Create Spike rasters of units by depth during brhavior state
if nargin<3
    sz = 15;
end

if flag
    disp('Analyzing run!')
    spikeTempL23 = squeeze(struct2cell(Spikes.spikes.L23Run)); % Spike temp for layer 23 of cell cell array
    spikeTempL23 = vertcat(spikeTempL23{:});
    spikeTempL5 = squeeze(struct2cell(Spikes.spikes.L5Run)); % Spike temp for layer 23 of cell cell array
    spikeTempL5 = vertcat(spikeTempL5{:});
else
    disp('Analyzing rest!')
    spikeTempL23 = squeeze(struct2cell(Spikes.spikes.L23Rest)); % Spike temp for layer 23 of cell cell array
    spikeTempL23 = vertcat(spikeTempL23{:});
    spikeTempL5 = squeeze(struct2cell(Spikes.spikes.L5Rest)); % Spike temp for layer 23 of cell cell array
    spikeTempL5 = vertcat(spikeTempL5{:});
end

% plot best representative neuron across trials
t = cellfun(@mean,spikeTempL23);
t(isnan(t)) = 0;
t = mean(t,1);
[~,idx] = max(t); %find max cell 
figure,hold on,ylim([0 50]),title('L2/3 single units')
trials = size(spikeTempL23,1);
rasterAvg = [];
rasterA = [];
Y = [];
for i = 1:30
    raster = spikeTempL23{i,idx};
    scatter(raster,i*ones(size(raster,1),1),sz,'filled','r')
    L23rasterAvg{i} = raster;
end
rasterAvg = vertcat(L23rasterAvg{:});
rasterAvg(isnan(rasterAvg)) = 1;
Y = discretize(rasterAvg,200); % Discretize into 10ms time bins by trial length (2000)/200
for i = 1:200
    rasterA(i) = length(rasterAvg(Y==i)); %Count number of spikes per bin
end
figure,plot(smoothdata(rasterA,'rloess',10),'r'),box off,ylim([0 25]),xlim([0 100]),title('L2/3 average response');
% spikeTemp = squeeze(struct2cell(Spikes.dspikes.L4Run)); % Spike temp for layer 23 of cell cell array
% spikeTemp = vertcat(spikeTemp{:});
% figure,hold on
% for i = 1:size(spikeTemp,1)
%     r = horzcat(spikeTemp(i,:));
%     raster{i} = vertcat(r{:});
%     plot(raster{i},(i/5)*ones(size(raster{i},1)),'.')
% end
% spikeTemp = [];


t = cellfun(@mean,spikeTempL5);
t(isnan(t)) = 0;
t = mean(t,1);
[~,idx] = max(t); %find max cell 
figure,hold on,ylim([0 50]),title('L5 single units')
rasterAvg = [];
raster = [];
Y = [];
rasterB = [];
for i = 1:trials
    raster = spikeTempL5{i,idx};
    scatter(raster,i*ones(size(raster,1),1),sz,'filled','b')
    L5rasterAvg{i} = raster;
end
rasterAvg = vertcat(L5rasterAvg{:});
rasterAvg(isnan(rasterAvg)) = 1;
Y = discretize(rasterAvg,200);
for i = 1:200
    rasterB(i) = length(rasterAvg(Y==i));
end
figure,plot(smoothdata(rasterB,'rloess',10),'b'),ylim([0 50]),xlim([0 100]), box off,title('L5 average response')
end
