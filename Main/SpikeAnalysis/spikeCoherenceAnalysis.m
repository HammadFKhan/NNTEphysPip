%% load in folders
foldername = strcat(uigetdir(pwd,'Input Directory'),'\');
filetype = 'mat'; % type of files to be processed
% Types currently supported .tif/.tiff, .h5/.hdf5, .raw, .avi, and .mat files
file = subdir(fullfile(foldername,['*.',filetype]));   % list of filenames (will search all subdirectories)
if isempty(file),disp('No Mat files where detected in this directory!'), return; end % handing for incorrect files

numFile = length(file);
%%
for fileNum = 1:numFile
    filename = file(fileNum).name;
    load(filename)
    %% Preps spike data for spike-spike coherence (since data is organized as spike x trial)
    fpath = pathname;
    path = [fpath,'/preAutoMerge'];
    Fs = 8192;
    % Read in kilosort data for matlab analysis
    SpikeClusters = readNPY(fullfile(path, 'spike_clusters.npy'));
    SpikeSamples = readNPY(fullfile(path, 'spike_times.npy'));
    %% Analysis
    Spikes.SpikeClusters = SpikeClusters; %Add one because of 0 index from python
    Spikes.SpikeSamples = SpikeSamples;
    Spikes = clusterSort(Spikes);
    Spikes = ISI(Spikes,0.01,Fs,0); %Spikes, Interval, Fs
    %%
    load chanMap % use for PFF
%     load UCLA_chanMap
    [spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] =...
        spikeTemplatePosition(fpath,ycoords);
    for i = 1:length(tempAmps)
        Spikes.Clusters(i).spikeDepth = templateDepths(i);
        Spikes.Clusters(i).spikeAmplitude = tempAmps(i);
        Spikes.Clusters(i).waveforms = waveforms(i,:);
        Spikes.Clusters(i).spikeDuration = templateDuration(i)/Fs*1000;
    end
    %% Select specific cells for coupling based on amplitude and waveform width
    spikeNum = [];
    for i = 1:length(Spikes.Clusters)
        spikeNum = [spikeNum;length(Spikes.Clusters(i).cluster)];
    end
    cells = find((templateDuration/Fs*1000)>1 & spikeNum>max(Spikes.SpikeSamples)/Fs/2);
    [val,I] = sort(tempAmps,'descend');
    for i = 1:length(I)
        if ismember(I(i),cells)
        else
            I(i) = nan;
        end
    end
    I(isnan(I)) = [];
    selectedCells = I(1:25); %Select top 10 best neurons
    for i = 1:length(Spikes.Clusters)
        if ismember(i,selectedCells)
            Spikes.Clusters(i).coherence = 1;
        else
            Spikes.Clusters(i).coherence = 0;
        end
    end
    C = nchoosek(selectedCells, 2);  % Indices combination for spikes
    %% Begin spike field coherence
    spikeCoherence = struct();
    spikea = {};
    spikeb = {};
    for spikeCom = 1:size(C,1)
        spikeCoherence(spikeCom).neuronA = C(spikeCom,1);
        spikeCoherence(spikeCom).neuronB = C(spikeCom,2);
        [trials,~] = makeSpikeWin(Spikes,C(spikeCom,1),Fs);
        spikeCoherence(spikeCom).spikea{1} = vertcat(trials{:});
        [trials,~] = makeSpikeWin(Spikes,C(spikeCom,2),Fs);
        spikeCoherence(spikeCom).spikeb{1} = vertcat(trials{:});
        disp(['Analyzing Spike ' num2str(C(spikeCom,1)) ' with Spike ' num2str(C(spikeCom,2))])
        pause(1)
        spikeCoherence(spikeCom).spikecoherence = ...
            tapered_spike_coherence_modified(spikeCoherence(spikeCom).spikea,spikeCoherence(spikeCom).spikeb,[],800,0);
    end
    %% Save data
    if ~exist([file(fileNum).folder '\output'],'dir')
        mkdir([file(fileNum).folder '\output']);
    end
    [folder_name,file_name,~] = fileparts(file(fileNum).name);
    if exist(fullfile([folder_name, '\output'],[file_name,'.mat']),'file')
        file_name = [file_name '_' datestr(now,30) '_'];
    end
    savepath = fullfile([folder_name, '\output'],[file_name,'.mat']);
    save(savepath, 'fpath','pathname','spikeCoherence','Spikes','Fs','-v7.3');
    clearvars -except file numFile fileNum filetype foldername
end
%% ----------- Plot data *optional* --------------- %%
plotspikeCoherence(spikeCoherence(1).spikea,spikeCoherence(1).spikeb,spikeCoherence(1).spikecoherence)


%% Some Functions
function [trials,win] = makeSpikeWin(Spikes,spikeId,Fs)
%% Begin analysis across clusters
winLen = 4;
win = double(0:winLen:max(Spikes.SpikeSamples)/Fs);
trials = cell(1,length(win));
% win(1:400) = [];
for i = 1:length(win)-1
    spiketm = Spikes.Clusters(spikeId).spikeTime-win(i);
    idx = find(spiketm>=0 & spiketm < winLen);
    trials{i} = zeros(1,length(0:0.001:winLen));
    spiked = discretize(spiketm(idx),0:0.001:winLen);
    if ~isnan(spiked)
        trials{i}(spiked) = 1;
    end
end
end

function plotspikeCoherence(spikea,spikeb,spikeCoherence)

interval = 0:2000;
spiker = spikea;

figure;
%****************
disp('Plotting rastergrams (slow) ...');
subplot('position',[0.1 0.4 0.4 0.55]);
make_nice_spike_raster(spiker);
V = axis;
axis([0 size(spikea{1},2) V(3) V(4)]);
grid on;
ylabel('Trial Number');
title(fprintf('Unit rasters'));
%*************
subplot('position',[0.1 0.1 0.4 0.3]);
smooth_window = 25;  % give sigma of 12.5ms
make_nice_mean_raster(spiker,smooth_window);
V = axis;
axis([0 size(spikea{1},2) V(3) V(4)]);
plot([interval(1),interval(1)],[V(3),V(4)],'k-'); hold on;
plot([interval(end),interval(end)],[V(3),V(4)],'k-'); hold on;
ylabel('Firing Rate');
xlabel('Time (ms)');

spiker = spikeb;
%****************
disp('Plotting rasters2 per trial (slow) ...');
subplot('position',[0.57 0.4 0.4 0.55]);
make_nice_spike_raster(spiker);
V = axis;
axis([0 size(spikeb{1},2) V(3) V(4)]);
grid on;
ylabel('Trial Number');
title(fprintf('Unit rasters'));
%*************
subplot('position',[0.57 0.1 0.4 0.3]);
smooth_window = 25;  % give sigma of 12.5ms
make_nice_mean_raster(spiker,smooth_window);
V = axis;
axis([0 size(spikeb{1},2) V(3) V(4)]);
plot([interval(1),interval(1)],[V(3),V(4)],'k-'); hold on;
plot([interval(end),interval(end)],[V(3),V(4)],'k-'); hold on;
ylabel('Firing Rate');
xlabel('Time (ms)');

figure;
colo = 'rbgy';
CNUM=1;
subplot(2,2,1);
for ii = 1:CNUM
    H = semilogx(spikeCoherence.freq{ii},spikeCoherence.coho{ii},[colo(ii),'-']); hold on;
    set(H,'Linewidth',2);
    H = semilogx(spikeCoherence.freq{ii},(spikeCoherence.coho{ii}+spikeCoherence.scoho{ii}),[colo(ii),':']); hold on;
    H = semilogx(spikeCoherence.freq{ii},(spikeCoherence.coho{ii}-spikeCoherence.scoho{ii}),[colo(ii),':']);
    H = semilogx(spikeCoherence.freq{ii},spikeCoherence.rcoho{ii},[colo(ii),'--']); hold on;
    set(H,'Linewidth',1);
end
ylabel('Coherence Magnitude');
xlabel('Frequency (Hz)');
title(sprintf('Spike A with Spike B'));

subplot(2,2,2);
for ii = 1:CNUM
    H = plot(spikeCoherence.freq{ii},spikeCoherence.phaso{ii},[colo(ii),'-']); hold on;
    set(H,'Linewidth',2);
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.phaso{ii}+spikeCoherence.sphaso{ii}),[colo(ii),':']); hold on;
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.phaso{ii}-spikeCoherence.sphaso{ii}),[colo(ii),':']);
    %H = plot(freq{ii},rphaso{ii},[colo(ii),'--']); hold on;
    set(H,'Linewidth',1);
end
ylabel('Coherence Angle');
xlabel('Freq');
title(sprintf('Spike A with Spike B'));

subplot(2,2,3);
for ii = 1:CNUM
    H = plot(spikeCoherence.freq{ii},spikeCoherence.spika_pow{ii},[colo(ii),'-']); hold on;
    set(H,'Linewidth',2);
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.spika_pow{ii}+spikeCoherence.sspika_pow{ii}),[colo(ii),':']); hold on;
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.spika_pow{ii}-spikeCoherence.sspika_pow{ii}),[colo(ii),':']);
end
ylabel('Spike Unit A Pow');
xlabel('Frequency (Hz)');
title(sprintf('Spike A with Spike B'));

subplot(2,2,4);
for ii = 1:CNUM
    H = plot(spikeCoherence.freq{ii},spikeCoherence.spike_pow{ii},[colo(ii),'-']); hold on;
    set(H,'Linewidth',2);
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.spike_pow{ii}+spikeCoherence.sspike_pow{ii}),[colo(ii),':']); hold on;
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.spike_pow{ii}-spikeCoherence.sspike_pow{ii}),[colo(ii),':']);
end
ylabel('Spike Unit B Pow');
xlabel('Frequency (Hz)');
title(sprintf('Spike B with Spike A'));

end
%****************** function to make a nice raster plot ****************
function make_nice_mean_raster(spmat,smooth_window)
%*********** spmat1 and spmat2 are spike matrices of two conditions you wish to compare
%*********** smooth_window ... gaussian smoothing in millisecs
numconds = size(spmat,2);
if (numconds==2)
    colo = [[1,0,0];[0,0,1]];
else
    colo = [[1,0,0];[0,0,1];[0,1,0];[1,1,0]];
end

for k = 1:numconds
    spud = spmat{1,k};
    numtrials(k) = size(spud,1);
    smorate = gauss_smooth(sum( spud(1:numtrials(k),:))/....
        numtrials(k),smooth_window)*1000;
    H = plot(smorate,'k-'); hold on;
    set(H,'Color',colo(k,:));
end

end


%******************* make a rastergram of the actual spikes on each trial
function make_nice_spike_raster(spmat)

colo = [[1,1,1];[0,0,0]];
colormap(colo);

totspike = [];
for cubo = 1:size(spmat,2)
    totspike = [totspike; (cubo*spmat{1,cubo}) ];
    totspike = [totspike; (cubo*spmat{1,cubo}) ];
    
end
totspike = totspike + 1;
imagesc(totspike);
V = axis;
axis([V(1) V(2) 0 size(totspike,1)]);

end

%**************************************************************
function output = gauss_smooth(input, window)
% Smoothing function:
% output = smooth(input, window)
% "Window" is the total kernel width.
% Input array must be one-dimensional.

input_dims = ndims(input);
input_size = size(input);
if input_dims > 2 | min(input_size) > 1,
    disp('Input array is too large.');
    return
end

if input_size(2) > input_size(1),
    input = input';
    toggle_dims = 1;
else
    toggle_dims = 0;
end

if window/2 ~= round(window/2),
    window = window + 1;
end
halfwin = window/2;

input_length = length(input);
%********* gauss window +/- 1 sigma
x = -halfwin:1:halfwin;
kernel = exp(-x.^2/(window/2)^2);
kernel = kernel/sum(kernel);

padded(halfwin+1:input_length+halfwin) = input;
padded(1:halfwin) = ones(halfwin, 1)*input(1);
padded(length(padded)+1:length(padded)+halfwin) = ones(halfwin, 1)*input(input_length);

output = conv(padded, kernel);
output = output(window:input_length+window-1);

if toggle_dims == 1,
    output = output';
end

end

