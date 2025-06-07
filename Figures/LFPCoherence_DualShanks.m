%% Dual shanks GP interareal LFP-LFP coherence
% Calculate LFP-LFP phase coherence across trials
if ~exist('LFP','var')
    load('Y:\Hammad\Ephys\LeverTask\DualShank\075356DualShank\Day3\Day3M2DualM1SingleRecording1_240727_180222\UCLA_chanmap_64F2\LFP.mat')
end

%%
dat1 = squeeze(LFP.probe1.genPhase.hitLFP{11}(1:4:end,:,:));
dat2 = squeeze(LFP.probe3.genPhase.hitLFP{11}(1:8:end,:,:));
figure,stack_plot(dat1,1,1,1000)
figure,imagesc(dat1)
figure,stack_plot(dat2,1,1,1000)
%% Plot it out
preCuephaseCorr = preCuephaseCorrHit;
postCuephaseCorr = postCuephaseCorrHit;

load myMap
f = figure;
subplot(321),imagesc(preCuephaseCorr),colormap(myMap),colorbar,
ylabel('Channel'),xlabel('Channel')
set(gca,'fontsize',14,'linewidth',1.5)
title('Pre-cue')%,caxis([0 1])
axis square

subplot(322),imagesc(postCuephaseCorr),colormap(myMap),colorbar,
ylabel('Channel'),xlabel('Channel')
set(gca,'fontsize',14,'linewidth',1.5)
title('Post-cue')%,caxis([0 1])
axis square

preCuephaseCorr = preCuephaseCorrMiss;
postCuephaseCorr = postCuephaseCorrMiss;

subplot(323),imagesc(preCuephaseCorr),colormap(myMap),colorbar,
ylabel('Channel'),xlabel('Channel')
set(gca,'fontsize',14,'linewidth',1.5)
title('Pre-cue')%,caxis([0 1])
axis square

subplot(324),imagesc(postCuephaseCorr),colormap(myMap),colorbar,
ylabel('Channel'),xlabel('Channel')
set(gca,'fontsize',14,'linewidth',1.5)
title('Post-cue')%,caxis([0 1])
axis square

preCuephaseCorr = preCuephaseCorrFA;
postCuephaseCorr = postCuephaseCorrFA;

subplot(325),imagesc(preCuephaseCorr),colormap(myMap),colorbar,
ylabel('Channel'),xlabel('Channel')
set(gca,'fontsize',14,'linewidth',1.5)
title('Pre-cue')%,caxis([0 1])
axis square

subplot(326),imagesc(postCuephaseCorr),colormap(myMap),colorbar,
ylabel('Channel'),xlabel('Channel')
set(gca,'fontsize',14,'linewidth',1.5)
title('Post-cue')%,caxis([0 1])
axis square
%% Do stats so we can pool across animals
preCuePCHit = mean(preCuephaseCorrHit);
postCuePCHit = mean(postCuephaseCorrHit);

preCuePCMiss = mean(preCuephaseCorrMiss);
postCuePCMiss = mean(postCuephaseCorrMiss);

preCuePCFA = mean(preCuephaseCorrFA);
postCuePCFA = mean(postCuephaseCorrFA);

%% Plot it out
figure,
subplot(131),customBoxplot([preCuePCHit',postCuePCHit'])
subplot(132),customBoxplot([preCuePCMiss',postCuePCMiss'])
subplot(133),customBoxplot([preCuePCFA',postCuePCFA'])
%% LOCAL FUNCTIONS
function [preCuephaseCorr, postCuephaseCorr] = getphaseLFP(xgp1,xgp2)
phaseLFP = cellfun(@(x) squeeze(x(:,:,1:1500)),xgp1,'UniformOutput',false);
preCuephaseLFP1 = angle(horzcat(phaseLFP{:}));
phaseLFP = cellfun(@(x) squeeze(x(:,:,1501:end)),xgp1,'UniformOutput',false);
postCuephaseLFP1 = angle(horzcat(phaseLFP{:}));

phaseLFP = cellfun(@(x) squeeze(x(:,:,1:1500)),xgp2,'UniformOutput',false);
preCuephaseLFP2 = angle(horzcat(phaseLFP{:}));
phaseLFP = cellfun(@(x) squeeze(x(:,:,1501:end)),xgp2,'UniformOutput',false);
postCuephaseLFP2 = angle(horzcat(phaseLFP{:}));

preCuephaseLFP = [preCuephaseLFP1;preCuephaseLFP2(1:2:end,:)];
postCuephaseLFP = [postCuephaseLFP1;postCuephaseLFP2(1:2:end,:)];
fprintf('Calculating lfp-lfp phase...\n')
NChan = size(preCuephaseLFP,1);
preCuephaseCorr = zeros(NChan,NChan);
postCuephaseCorr = zeros(NChan,NChan);
for i = 1:NChan
    for j = 1:NChan
        preCuephaseCorr(i,j) = circ_corrcc(preCuephaseLFP(i,:),preCuephaseLFP(j,:));
        postCuephaseCorr(i,j) = circ_corrcc(postCuephaseLFP(i,:),postCuephaseLFP(j,:));
    end
    disp(['Chan: ' num2str(i)])
end
fprintf('done\n')
end