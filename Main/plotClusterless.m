function plotClusterless(Spikes,Ripples,filtData,Intan)
amplifier_data = Intan.amplifier_data;

set(0,'DefaultFigureWindowStyle','docked')
addpath(genpath('Figures'));
disp('Plotting...')

try 
figure,spikePlot = Show_Spikes(Spikes.binary);
catch ME
end

try 
figure('Name', 'Unfiltered Data'),stack_plot(amplifier_data);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
    print(gcf,'-painters','-depsc', 'Figures/Raw_data.eps', '-r250');
catch ME
end

try
figure('Name','Singe Unit Waveforms'),SingleUnits(Spikes.allSpikes);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
    print(gcf,'-painters','-depsc', 'Figures/Single_unit_waveforms.eps', '-r250');
catch ME
end

try
figure('Name','Multi-Unit Activity'),MUA(filtData.MUA);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
    print(gcf,'-painters','-depsc', 'Figures/MultiUnit.eps', '-r250');
catch ME
end

try
figure('Name','LFP'),LFP(filtData.LFP); 
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
    print(gcf,'-painters','-depsc', 'Figures/LFP.eps', '-r250');
catch ME
end

try 
ripplePlot(Ripples);
figure('Name','Ripple Stimulus'),h = htmp(correlation-mean(correlation,'all'),50);caxis([0 1]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
    print(gcf,'-painters','-depsc', 'Figures/Peristimulus.eps', '-r250');
catch ME
end