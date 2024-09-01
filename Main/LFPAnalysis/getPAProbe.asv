function [PA] = getPAProbe(LFP,IntanBehaviour,z_score,nPerm,plotFlag,parameters)

% Reference : Spontaneous travelling cortical waves gate perception in 
% behaving primates, Nature, 2020 

% Similair to getPA but adjusted for silicon probe data. (HK 06/2024)
parameters.rows = size(LFP.hitxgp{1},1);
parameters.cols = size(LFP.hitxgp{1},2);
nElectrodes = parameters.rows*parameters.cols;

xgpHit = LFP.hitxgp;
xgpMiss = LFP.missxgp;
xgpMIhit = LFP.MIhitxgp;
xgpMIFA = LFP.MIFAxgp;

[PA.Hit,PA.Hit_ang] = calPhaseAlignment(xgpHit);
[PA.Miss,PA.Miss_ang] = calPhaseAlignment(xgpMiss);
[PA.MIHit,PA.MIHit_ang] = calPhaseAlignment(xgpMIhit);
[PA.MIFA,PA.MIFA_ang] = calPhaseAlignment(xgpMIFA);

if plotFlag == 1
    figure();
    title("Phase Alignment across all electrodes - Hits")
    subplot(2,1,1);
    imagesc(IntanBehaviour.cueHitTrace(1).time,1:nElectrodes,reshape(PA.Hit,[],size(PA.Hit,3))); colormap(hot);
    ylabel("Electrodes");xlabel("Time (s)");
    xline(0,'-w','Cue','LabelVerticalAlignment','top');
    subplot(2,1,2);
    title("Phase Alignment across all electrodes - Misses")
    imagesc(IntanBehaviour.cueMissTrace(1).time,1:nElectrodes,reshape(PA.Miss,[],size(PA.Miss,3))); colormap(hot);
    ylabel("Electrodes");xlabel("Time (s)");
    xline(0,'-w','Cue','LabelVerticalAlignment','top');
    % subplot(3,1,3);
    % title("Phase Alignment across all electrodes - False Alarms")
    % imagesc(IntanBehaviour.cueHitTrace(1).time,1:32,reshape(PAFA,[],size(PAFA,3))); colormap(hot);
    % ylabel("Electrodes");xlabel("Time (s)");
    % xline(0,'-w','Cue','LabelVerticalAlignment','top');
    
    figure();
    title("Phase Alignment averaged across Electrodes")
    plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PA.Hit,[1 2])),'-r','LineWidth',1.2); hold on;
    plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PA.Miss,[1 2])),'-k','LineWidth',1);
    % plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PAFA,[1 2])),'-k','LineWidth',1);
    ylabel("Phase Alignment"); xlabel("Time (s)");
    xline(0,'--r','Cue','LabelVerticalAlignment','top');
    xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
    title('Phase Alignment for Hits');box off;legend('Hits','Miss');
end

% z-scoring
if z_score == 1
    xgpComb = [xgpHit xgpMiss];
    nHit = size(xgpHit,2);
    nMiss = size(xgpMiss,2);
    nTot = nHit + nMiss;
    
    xgpComb2 = [xgpMIhit xgpMIFA];
    nMIhit = size(xgpMIhit,2);
    nMIFA = size(xgpMIFA,2);
    nTot2 = nMIhit + nMIFA;
    
    nullDistHit = zeros(parameters.rows,parameters.cols,size(PA.Hit,3),nPerm);
    nullDistMiss = zeros(parameters.rows,parameters.cols,size(PA.Miss,3),nPerm);
    
    nullDistMIhit = zeros(parameters.rows,parameters.cols,size(PA.MIHit,3),nPerm);
    nullDistMIFA= zeros(parameters.rows,parameters.cols,size(PA.MIFA,3),nPerm);
    
    tic
    parfor j=1:nPerm
        randIndex = randperm(nTot);
        xgpHitRand = xgpComb(randIndex(1:nHit));
        xgpMissRand = xgpComb(randIndex(nHit+1:end));
        nullDistHit(:,:,:,j) = calPhaseAlignment(xgpHitRand);
        nullDistMiss(:,:,:,j) = calPhaseAlignment(xgpMissRand);
        
        randIndex2 = randperm(nTot2);
        xgpMIhitRand = xgpComb2(randIndex2(1:nMIhit));
        xgpMIFARand = xgpComb2(randIndex2(nMIhit+1:end));
        nullDistMIhit(:,:,:,j) = calPhaseAlignment(xgpMIhitRand);
        nullDistMIFA(:,:,:,j) = calPhaseAlignment(xgpMIFARand);
        j
    end
    toc
    PA.muHit = mean(nullDistHit,4); % Mean of the null distribution
    PA.sigmaHit = std(nullDistHit,0,4); % Standard deviation of null distribution
    PA.muMiss = mean(nullDistMiss,4); % Mean of the null distribution
    PA.sigmaMiss = std(nullDistMiss,0,4); % Standard deviation of null distribution
    
    PA.muMIhit = mean(nullDistMIhit,4); % Mean of the null distribution
    PA.sigmaMIhit  = std(nullDistMIhit,0,4); % Standard deviation of null distribution
    PA.muMIFA = mean(nullDistMIFA,4); % Mean of the null distribution
    PA.sigmaMIFA = std(nullDistMIFA,0,4); % Standard deviation of null distribution
    
    PA.Hitz = (PA.Hit-PA.muHit)./PA.sigmaHit;
    PA.Missz = (PA.Miss-PA.muMiss)./PA.sigmaMiss;
    
    PA.MIhitz = (PA.MIHit-PA.MIHit)./PA.sigmaMIhit;
    PA.MIFAz = (PA.MIFA-PA.muMIFA)./PA.sigmaMIFA;
    
    if plotFlag == 1
        figure();
        plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PA.Hitz,[1 2])),'-r','LineWidth',1.2); hold on;
        plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PA.Missz,[1 2])),'-k','LineWidth',1);
        ylabel("z-score"); xlabel("Time (s)");
        xline(0,'--r','Cue','LabelVerticalAlignment','top');
        xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
        title('z-scored Phase Alignment averaged across electrodes');box off; legend('Hits','Miss');  
        
        figure();
        subplot(2,1,1);
        title("Phase Alignment across all electrodes - Hits")
        imagesc(IntanBehaviour.cueHitTrace(1).time,1:nElectrodes,reshape(PA.Hitz,[],size(PA.Hitz,3))); colormap(hot);
        ylabel("Electrodes");xlabel("Time (s)");
        xline(0,'-w','Cue','LabelVerticalAlignment','top');caxis([-4 8]);
        xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
        subplot(2,1,2);
        title("Phase Alignment across all electrodes - Misses")
        imagesc(IntanBehaviour.cueMissTrace(1).time,1:nElectrodes,reshape(PA.Missz,[],size(PA.Missz,3))); colormap(hot);
        ylabel("Electrodes");xlabel("Time (s)");
        xline(0,'-w','Cue','LabelVerticalAlignment','top');caxis([-4 8]);
        yline(2.56);
    end
end