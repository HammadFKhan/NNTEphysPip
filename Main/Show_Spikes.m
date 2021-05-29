function spikePlot  = Show_Spikes(Spikes)
y = length(Spikes(:,1));
x = length(Spikes(1,:));
spikePlot = zeros(y,x);
baseline = 0;
for i = 1:y
    for ii = 1:x
        if Spikes(i,ii) == 1
            spikePlot(i,ii) = Spikes(i,ii)+baseline;
        else
            spikePlot(i,ii) = NaN;
        end
    end
    plot(spikePlot(i,:),'.','MarkerEdgeColor','k'); hold on;
%     plot_raster(spikePlot(i,:)); hold on;
    baseline = baseline + 1;
end
ylim ([0,y+1]);
end
