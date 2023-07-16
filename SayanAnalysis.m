addpath(genpath('main'));
foldername = strcat(uigetdir(pwd,'Input Directory'),'\');
filetype = 'mat'; % type of files to be processed
% Types currently supported .tif/.tiff, .h5/.hdf5, .raw, .avi, and .mat files
file = subdir(fullfile(foldername,['*.',filetype]));   % list of filenames (will search all subdirectories)
if isempty(file),disp('No Mat files where detected in this directory!'), return; end % handing for incorrect files
numFile = length(file);
for fileNum = 1:numFile
    filename = file(fileNum).name;
    load(filename)
    Spikes = spikeRateAnalysis(Spikes);
    L23SR{fileNum} = Spikes.spikeRate.mL23SR;
    L5SR{fileNum} = Spikes.spikeRate.mL5SR;
end
%%
[dN,t]=binspikes(Spikes.Clusters(3).spikeTime,100);
params.Fs = 8192;
params.tapers = [5 9];
params.fpass = [1 80];
[C,phi,S12,S1,S2,f]=coherencypt(Spikes.Clusters(3).spikeTime,Spikes.Clusters(4).spikeTime(1:16750),params);
%%
binFR = [];
neuronId = 122;
figure,
histogram(Spikes.Clusters(neuronId).ISI,0:.01:2)
figure,
for i = 200:250
    spike = Spikes.Clusters(neuronId).spikeTime-Spikes.Clusters(neuronId).spikeTime(1);
    sTemp = spike-(i-1)*2;
    idx = find(sTemp>0 & sTemp<=2);
    blah = discretize(sTemp(idx),0:.01:2)';
    binFR = vertcat(binFR,blah);
    plot(sTemp(idx),i*ones(size(idx)),'b.'),hold on
end