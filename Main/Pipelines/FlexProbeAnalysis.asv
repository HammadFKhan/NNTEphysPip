%% Flex probe pipeline
% Alterations in kilosort layout and changes in spike sorting metrics

addpath(genpath('chronux'));
addpath(genpath('Kilosort'));
addpath(genpath('npy-matlab'));
addpath(genpath('spikes-master'));
IntanConcatenateFlex
savepath = fullfile(pathname,['loadme','.mat']);
save(savepath,'ds_filename');
clearvars -except ds_filename
%%
% Intan = read_Intan_RHD2000_file(); %load intan data
% ParpoolConfig
fpath    = Intan.path; % where on disk do you want the analysis? ideally and SSD...
pathToYourConfigFile = strcat(pwd,'/main/'); % for this example it's ok to leave this path inside the repo, but for your own config file you *must* put it somewhere else!  
run(fullfile(pathToYourConfigFile, 'config_FlexProbe.m'))
make_Flex4BSqChannelMap(fpath,s); % Creates channel map for electrode array
% make_UCLAMouseChannelMap(fpath);
kilosortPrep(Intan.allIntan,fpath);
set(0,'DefaultFigureWindowStyle','normal')
rez = KilosortAnalysis(fpath,ops);
% now fire up Phy and check these results. There should still be manual
% work to be done (mostly merges, some refinements of contaminated clusters). 
%% AUTO MERGES 
% after spending quite some time with Phy checking on the results and understanding the merge and split functions, 
% come back here and run Kilosort's automated merging strategy. This block
% will overwrite the previous results and python files. Load the results in
% Phy again: there should be no merges left to do (with the default simulation), but perhaps a few splits
% / cleanup. On realistic data (i.e. not this simulation) there will be drift also, which will usually
% mean there are merges left to do even after this step. 
% Kilosort's AUTO merges should not be confused with the "best" merges done inside the
% benchmark (those are using the real ground truth!!!)
data = matfile(ds_filename);
% LFP
set(0,'DefaultFigureWindowStyle','normal')
% LFP = fastpreprocess_filtering(flip(Intan.allIntan,1),8192); %Only run for PFF data
temp = data.Intan;
LFP = fastpreprocess_filtering(temp.allIntan,2000);
LFP = bestLFP(LFP);
LFP = bandFilter(LFP,'depth'); % Extract LFPs based on 'depth' or 'single'
%%
set(0,'DefaultFigureWindowStyle','normal')
% load chanMap % use for PFF
load Flex4BSq_Chan_Map
[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] =...
    spikeTemplatePosition(fpath,ycoords);
waveformsFixed = [];
for i = 1:size(waveforms,1)
    t = isnan(waveforms(i,:));
    t = find(t==1);
    temp = waveforms(i,:);
    temp(t) = 0;
    x1 = interp1(1:length(temp),temp,1:0.5:length(temp),'pchip');
    waveformsFixed(i,:) = x1;
end
figure,plot(waveformsFixed')
[X] = featureProject(waveformsFixed',1,0);
