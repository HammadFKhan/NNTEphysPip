function [Ripples,filtData] = SWR(Intan)
%% Intan Data
% Parse
amplifier_data = Intan.amplifier_data;
t_amplifier = Intan.t_amplifier;
addpath(genpath('main'));
% load chanMap.mat
% for i = 1:length(chanMap)
%     amplifier_data(i,:) = amplifier_data(chanMap(i),:);
% end
%% Data Processing
filtData = preprocess_filtering(amplifier_data,t_amplifier);
try
    Ripples = rippleDetection(filtData);
catch ME
    disp('Error in Ripple detection!')
end

try 
    Ripples = rippleAnalysis(filtData,Ripples);
catch ME
    disp('Error in Ripple Analysis!')
end
%%
% data = filtData.lowpassData';
% spacing = 2E-5;
% [CSDoutput]  = CSD(data,20000,spacing,'inverse',spacing*5);

