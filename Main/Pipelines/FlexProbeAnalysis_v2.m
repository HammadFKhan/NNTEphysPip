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
% LFP = bestLFP(LFP);
% LFP = bandFilter(LFP,'depth'); % Extract LFPs based on 'depth' or 'single'
%%
addpath(genpath('C:\Users\khan332\Documents\GitHub\generalized-phase'));
xo = bandpass_filter(LFP.LFP,5,40,1000); %x,f1,f2,Fs
sz = size(xo);
xo = reshape(xo,sz(1),1,sz(2));
%% Feed in chunks of data for GP so that it actually runs decently
sz = size(xo);
chunksz = floor(sz(3)/100);
for n = 1:100
    win = (n-1)*chunksz+1:n*chunksz;
    xgp{n} = generalized_phase(xo(:,:,win),1000, 0 ,0);
    n
end
xgp{n+1} = generalized_phase(xo(:,:,win(end):end),1000, 0 ,1);
%%
NChan = sz(1);
for i = 1:NChan
    phaseLFP = cellfun(@(x) x(i,:), xgp, 'UniformOutput', false);
    phaseLFP = angle(horzcat(phaseLFP{:}));
end
fprintf('done\n')

fprintf('Calculating lfp-lfp phase...')
phaseCorr = zeros(NChan,NChan);
phaseLFP = cellfun(@squeeze,xgp,'UniformOutput',false);
phaseLFP = angle(horzcat(phaseLFP{:}));
for i = 1:NChan
    for j = 1:NChan
        phaseCorr(i,j) = circ_corrcc(phaseLFP(i,1:100000),phaseLFP(j,1:100000));
    end
end
fprintf('done\n')