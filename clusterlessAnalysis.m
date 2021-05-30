function [Spikes,Ripples] = clusterlessAnalysis(Intan)
%% Intan Data
% Parse
amplifier_data = Intan.amplifier_data;
t_amplifier = Intan.t_amplifier;
path = Intan.path;
filename = Intan.filename;
kilosortPrep(amplifier_data,path,filename) %Prep raw data for Kilosort Analysis
addpath(genpath('main'));
load chanMap.mat
for i = 1:length(chanMap)
    amplifier_data(i,:) = amplifier_data(chanMap(i),:);
end
%% Data Processing
filtData = preprocess_filtering(amplifier_data,t_amplifier);
Spikes = spikeSorting(filtData);
Ripples = rippleDetection(filtData);
Ripples = rippleAnalysis(filtData,Ripples,Spikes);
%%
PeriStimt = sum(Ripples.rippleOnset.PeriStim,3);
[vectorized,~] = cosine_similarity(PeriStimt(:,1:5000),50);
correlation = abs(corr(vectorized));
correlation(isnan(correlation)) = 0;
%%
% data = filtData.lowpassData';
% spacing = 2E-5;
% [CSDoutput]  = CSD(data,20000,spacing,'inverse',spacing*5);

%% Plots
set(0,'DefaultFigureWindowStyle','docked')
addpath(genpath('Figures'));
disp('Plotting...')
figure,spikePlot = Show_Spikes(Spikes.binary);
figure('Name', 'Unfiltered Data'),stack_plot(amplifier_data);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
    print(gcf,'-painters','-depsc', 'Figures/Raw_data.eps', '-r250');
figure('Name','Singe Unit Waveforms'),SingleUnits(Spikes.allSpikes);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
    print(gcf,'-painters','-depsc', 'Figures/Single_unit_waveforms.eps', '-r250');
figure('Name','Multi-Unit Activity'),MUA(filtData.MUA);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
    print(gcf,'-painters','-depsc', 'Figures/MultiUnit.eps', '-r250');
figure('Name','LFP'),LFP(filtData.LFP); 
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
    print(gcf,'-painters','-depsc', 'Figures/LFP.eps', '-r250');
ripplePlot(Ripples);
figure('Name','Ripple Stimulus'),h = htmp(correlation-mean(correlation,'all'),50);caxis([0 1]);