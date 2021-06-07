function [Spikes,Ripples,filtData] = clusterlessAnalysis(Intan)
%% Intan Data
% Parse
amplifier_data = Intan.amplifier_data;
t_amplifier = Intan.t_amplifier;
path = Intan.path;
filename = Intan.filename;
kilosortPrep(amplifier_data,path); %Prep raw data for Kilosort Analysis
addpath(genpath('main'));
% load chanMap.mat
% for i = 1:length(chanMap)
%     amplifier_data(i,:) = amplifier_data(chanMap(i),:);
% end
%% Data Processing
filtData = preprocess_filtering(amplifier_data,t_amplifier);
Spikes = spikeSorting(filtData);
try
    Ripples = rippleDetection(filtData);
catch ME
    disp('Error in Ripple detection!')
end

try 
    Ripples = rippleAnalysis(filtData,Ripples,Spikes);
catch ME
    disp('Error in Ripple Analysis!')
end
%%
% data = filtData.lowpassData';
% spacing = 2E-5;
% [CSDoutput]  = CSD(data,20000,spacing,'inverse',spacing*5);

%% Plots
plotClusterless(Spikes,Ripples,filtData,Intan)