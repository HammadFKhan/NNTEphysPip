function [SqI, PE, TS] = SeqIndexDB(Data,NumEntropyBins)
%The Sequentiality Index (SqI) is is the product of two entropy measures:
% entropy of the distribution of peak activity time (Peak Entropy) and entropy of the
% spatsity of neural activity (Temporal Sparsity). The peak entropy will be maximal if the
% there is a uniform distribution of neurons whose peaks tile the entire
% delay period. One minus the sparsiity entropy will me maxminal when the activity at
% each moment comes from a single neuron. Thus SqI is maximal when the
% neural sequence is a perfect diagonal line across the entire window of
% interest with no overalap between neurons (a perfect "neural chain")
% Inputs:
% Data: a 3D matrix of neural activity in which index of trials, neurons,
% and time are in 1, 2, 3 dimension repectively
% NumEntropyBins: the number of entrpopy bins used to estimate the
% distribution of peak times and the Peak Entropy.

% Outputs:
% SqI: Sequentiality Index
% PE: Peak Entropy
% TS: mean Temporal Sparsity
% Dean Buonomano and Shanglin Zhou 9/14/2019

[numTrials,numUnits, max_t]=size(Data);
binsize = (max_t/(NumEntropyBins-1));


for i = 1:numTrials
    data_temp = squeeze(Data(i,:,:));% data for single trial
    
    % PEAK ENTROPY
    [peakval, peaktimes]=max(data_temp,[],2);
    [N, X]=hist(peaktimes,0:binsize:max_t);
    probpeak = N./numUnits;
    probpeak = probpeak(probpeak>0);
    pe = sum(-probpeak.*log2(probpeak));
    PE(i) = pe/log2(NumEntropyBins);  %Normalized Peak Entropy
    
    % MEAN TEMPORAL SPARSITY
    dn = data_temp./sum(data_temp); % normalized activity at each time
    ts = nansum(-dn.*log2(dn));
    ts = mean(ts);
    TS(i) = 1 - ts./log2(numUnits); %normalized temporal sparsity
    
    SqI(i) = sqrt(PE(i)*TS(i));
end
SqI = mean(SqI);
PE = mean(PE);
TS = mean(TS);



