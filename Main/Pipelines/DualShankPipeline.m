%% Duals shank recording analysis pipeline. 
% Here we are just porting in a bunch of useful function from the original
% pipeline but making it exclusive for dual shank recordings. Want to do
% this for the sake of time and brevity.
%% Custom call of loading in dual shank data with active electrodes. 
% Here we need to make sure that the electrodes we map are correct based on
% the file that is being imported. 
addpath(genpath('Main'));
intandsFlag = 1; %make LFPs
activeElectrodes = 1:64;
%chanMapFile = 'UCLA_chanmap_fixed.mat'; %UCLA Sharp
chanMapFile = 'UCLA_chanmap_64F2.mat';
if ~exist('pathname','var')
    pathname = uigetdir(pwd,'Input Directory');
end
ds_filename1 = intanPreprocessingDualShanks(pathname,chanMapFile,intandsFlag,activeElectrodes);
activeElectrodes = 65:128;
chanMapFile = 'UCLA_chanmap_fixed.mat'; %UCLA Sharp
%chanMapFile = 'UCLA_chanmap_64F2.mat';
ds_filename2 = intanPreprocessingDualShanks(pathname,chanMapFile,intandsFlag,activeElectrodes);

