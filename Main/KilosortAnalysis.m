function rez = KilosortAnalysis(fpath,ops)
%% Kilosort Analysis

load chanMap.mat %### load in matlab file for channel configuration



% This part runs the normal Kilosort processing on the simulated data
[rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)
% This runs the benchmark script. It will report both 1) results for the
% clusters as provided by Kilosort (pre-merge), and 2) results after doing the best
% possible merges (post-merge). This last step is supposed to
% mimic what a user would do in Phy, and is the best achievable score
% without doing splits. 
benchmark_simulation(rez, fullfile(pwd, 'eMouseGroundTruth.mat'));

% save python results file for Phy
mkdir(fpath,'preAutoMerge')
rezToPhy(rez, [fpath,'/preAutoMerge']);

fprintf('Kilosort took %2.2f seconds \n', toc)

rez2 = merge_posthoc2(rez);
disp('Automerging completed!')
benchmark_simulation(rez2, fullfile(pwd, 'eMouseGroundTruth.mat'));

% save python results file for Phy
mkdir(fpath,'postAutoMerge')
rezToPhy(rez2, [fpath,'/postAutoMerge']);

% remove temporary file
delete(ops.fproc);

