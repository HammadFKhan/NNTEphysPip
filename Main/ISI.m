function Spikes = ISI(Spikes,Interval,Fs)
if nargin<3 || strcmp(Fs,''), Fs = 20000; end
if nargin<2 || strcmp(Interval,''), Interval = 0.005;disp(['Interval set as ' Interval 'ms']); end
sizePlot = ceil(sqrt(length(Spikes.Clusters)));
figure('name','ISI')
for i = 1:length(Spikes.Clusters)
    cluster = double(Spikes.Clusters(i).cluster)/Fs;
    ISI = diff(cluster);
    Spikes.Clusters(i).ISI = ISI;
    Spikes.Clusters(i).spikeTime = cluster;
    subplot(sizePlot,sizePlot,i),histogram(Spikes.Clusters(i).ISI,'BinWidth',Interval,'FaceColor',[0 0 0]),axis tight,box off;
    title(['Cluster ' num2str(i)]);
%     set(gca,'YTick','','YTickLabel','');
    xlabel('Time (s)')
    xlim([0 1]);
end

    