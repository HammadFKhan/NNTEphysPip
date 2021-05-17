%% Intan Data
clear; clc; 
close all;
read_Intan_RHD2000_file
% % amplifier_data = amplifier_data(:,1:100000);
% t_amplifier = t_amplifier(1:100000);
amplifier_data_sorted = channelSortEdge(amplifier_data);
amplifier_data = amplifier_data_sorted;
%%
% set(0,'DefaultFigureWindowStyle','docked')
filtData = preprocess_filtering(amplifier_data,t_amplifier);
Spikes = spikeSorting(filtData);
Ripples = rippleDetection(filtData);
Ripples = rippleAnalysis(filtData,Ripples);
ripplePlot(Ripples);
%% Plots
disp('Plotting...')
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
%%
try
encoder_data = convert_encoder(board_adc_data(2,:),timestamps);
catch ME
    disp('No ADC data')
    return
end
figure('Name','Pulse Data');plot(encoder_data.rotate_pulse);
figure('Name','Angular Distance');bar(encoder_data.ang_distance);
figure('Name','Angular Velocity');bar(encoder_data.ang_velocity,'FaceColor',[.16 .835 .384],'EdgeColor','none');
figure('Name','Avg. Angular Velocity');avgV = movmean(encoder_data.ang_velocity,2);bar(avgV,'FaceColor',[.16 .835 .384],'EdgeColor','none');
