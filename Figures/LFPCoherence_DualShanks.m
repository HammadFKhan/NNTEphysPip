%% Dual shanks GP interareal LFP-LFP coherence
% Calculate LFP-LFP phase coherence across trials
if ~exist('LFP','var')
    load('Y:\Hammad\Ephys\LeverTask\DualShank\075356DualShank\Day3\Day3M2DualM1SingleRecording1_240727_180222\UCLA_chanmap_64F2\LFP.mat')
end
xgp1 = LFP.probe1.genPhase.hitxgp;
xgp2 = LFP.probe2.genPhase.hitxgp;
xgp3 = LFP.probe3.genPhase.hitxgp;

[preCuephaseCorrHit, postCuephaseCorrHit] = getphaseLFP(xgp1,xgp3);
xgp1 = LFP.probe1.genPhase.missxgp;
xgp3 = LFP.probe3.genPhase.missxgp;
[preCuephaseCorrMiss, postCuephaseCorrMiss] = getphaseLFP(xgp1,xgp3);

xgp1 = LFP.probe1.genPhase.MIFAxgp;
xgp3 = LFP.probe3.genPhase.MIFAxgp;
[preCuephaseCorrFA, postCuephaseCorrFA] = getphaseLFP(xgp1,xgp3);
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
preCuePCHit = mean(mean(preCuephaseCorrHit));
postCuePCHit = mean(mean(postCuephaseCorrHit));

preCuePCMiss = mean(mean(preCuephaseCorrMiss));
postCuePCMiss = mean(mean(postCuephaseCorrMiss));

preCuePCFA = mean(mean(preCuephaseCorrFA));
postCuePCFA = mean(mean(postCuephaseCorrFA));

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