function neuralTrajAnalysis(Spikes,Waves,Behaviour)
%% Take orthoganal latent dimensions and do statistics across trials
X = arrayfun(@(x) vertcat(x.xorth),Spikes.GPFA.seqTrainHit,'UniformOutput',false);
X = horzcat(X{:});
neuralTrajHit = reshape(X,size(X,1),Spikes.GPFA.seqTrainHit(1).T,[]);

X = arrayfun(@(x) vertcat(x.xorth),Spikes.GPFA.seqTrainMiss,'UniformOutput',false);
X = horzcat(X{:});
neuralTrajMiss = reshape(X,size(X,1),Spikes.GPFA.seqTrainMiss(1).T,[]);

X = arrayfun(@(x) vertcat(x.xorth),Spikes.GPFA.seqTrainMIHit,'UniformOutput',false);
X = horzcat(X{:});
neuralTrajMIHit = reshape(X,size(X,1),Spikes.GPFA.seqTrainMIHit(1).T,[]);

X = arrayfun(@(x) vertcat(x.xorth),Spikes.GPFA.seqTrainMIFA,'UniformOutput',false);
X = horzcat(X{:});
neuralTrajMIFA = reshape(X,size(X,1),Spikes.GPFA.seqTrainMIFA(1).T,[]);

X = arrayfun(@(x) vertcat(x.xorth),Spikes.GPFA.seqTrainHitMiss,'UniformOutput',false);
X = horzcat(X{:});
neuralTrajHitMiss = reshape(X,size(X,1),Spikes.GPFA.seqTrainHitMiss(1).T,[]);

X = arrayfun(@(x) vertcat(x.xorth),Spikes.GPFA.seqTrainMIHitFA,'UniformOutput',false);
X = horzcat(X{:});
neuralTrajMIHitFA = reshape(X,size(X,1),Spikes.GPFA.seqTrainMIHitFA(1).T,[]);

%% Calculate average trajectories and divergence based on trial difference
% local function call for meaning based on combined PCA of trial conditions
X = neuralTrajHitMiss;
hittrials = 1:length(Behaviour.cueHitTrace);
misstrials = length(Behaviour.cueHitTrace)+1:size(X,3);
rh = meanTraj(X,hittrials,6)'; %trajectory variable and predefined conditional trial indexes
rm = meanTraj(X,misstrials,6)'; %trajectory variable and predefined conditional trial indexes
% rhminit = abs(rh(1,:)-rm(1,:));
% rm = rm-rhminit;
t = linspace(-Behaviour.parameters.windowBeforeCue,Behaviour.parameters.windowAfterCue,size(rh,1));
% Find important indices in array
stimStart = interp1(t,1:length(t),0,'nearest'); % Find zero of data which reflects some task condition we made
rawreactionTime = arrayfun(@(x) x.reactionTime,Behaviour.cueHitTrace);
mreactionTime = interp1(t,1:length(t),mean(rawreactionTime),'nearest');
reactionTime = interp1(t,1:length(t),rawreactionTime,'nearest');
reactionTime(isnan(reactionTime)) = length(t);
if exist('Waves','var')
    rawWaveDensityhit = arrayfun(@(x) vertcat(x.wavePresent),Waves.wavesHit,'UniformOutput',false);rawWaveDensityhit= vertcat(rawWaveDensityhit{:});
    rawWavePGDhit = arrayfun(@(x) vertcat(x.PGD), Waves.wavesHit,'UniformOutput',false);rawWavePGDhit = vertcat(rawWavePGDhit{:});
    rawWaveSpeedhit = arrayfun(@(x) vertcat(x.s), Waves.wavesHit,'UniformOutput',false);rawWaveSpeedhit = vertcat(rawWaveSpeedhit{:});
    
    rawWaveDensitymiss = arrayfun(@(x) vertcat(x.wavePresent),Waves.wavesMiss,'UniformOutput',false);rawWaveDensitymiss= vertcat(rawWaveDensitymiss{:});
    rawWavePGDmiss = arrayfun(@(x) vertcat(x.PGD), Waves.wavesMiss,'UniformOutput',false);rawWavePGDmiss = vertcat(rawWavePGDmiss{:});
    rawWaveSpeedmiss = arrayfun(@(x) vertcat(x.s), Waves.wavesMiss,'UniformOutput',false);rawWaveSpeedmiss = vertcat(rawWaveSpeedmiss{:});
    
    waveDensityhit = [];
    wavePGDhit = [];
    waveSpeedhit = [];
    
    waveDensitymiss = [];
    wavePGDmiss = [];
    waveSpeedmiss = [];
    
    win = ceil(1:20:size(rawWaveDensityhit,2));
    for n = 1:length(win)-1
        waveDensityhit = horzcat(waveDensityhit,sum(rawWaveDensityhit(:,win(n):win(n+1)),2)); %convert to density/sec
        wavePGDhit = horzcat(wavePGDhit,mean(rawWavePGDhit(:,win(n):win(n+1)),2));
        waveSpeedhit = horzcat(waveSpeedhit,mean(rawWaveSpeedhit(:,win(n):win(n+1)),2));
        waveDensitymiss = horzcat(waveDensitymiss,sum(rawWaveDensitymiss(:,win(n):win(n+1)),2)); %convert to density/sec
        wavePGDmiss = horzcat(wavePGDmiss,mean(rawWavePGDmiss(:,win(n):win(n+1)),2));
        waveSpeedmiss = horzcat(waveSpeedmiss,mean(rawWaveSpeedmiss(:,win(n):win(n+1)),2));
    end
    
end
% waveDensity = interp1(t,1:length(t),mean(rawWaveDensity,1),'nearest');
rawrewardTime = cell2mat(arrayfun(@(x) x.rewardIndex - x.LFPIndex(1501),Behaviour.cueHitTrace,'UniformOutput',false))/1000;
mrewardTime = interp1(t,1:length(t),mean(rawrewardTime),'nearest');
mrewardTime(isnan(reactionTime)) = length(t);

% Calculate differences in trajectories r'c(t)/||r'c(t)||
[neuralsimhitmiss,neuraldiffhitmiss,rprimeh,rprimem] = neuralTrajDiff(rh,rm);

%%% Do analysis as a function of reaction time
threshRT = median(rawreactionTime);
slowTrials = rawreactionTime>threshRT;
fastTrials = rawreactionTime<threshRT;
mreactionTimeslow = interp1(t,1:length(t),mean(rawreactionTime(slowTrials)),'nearest');
mreactionTimefast = interp1(t,1:length(t),mean(rawreactionTime(fastTrials)),'nearest');
rf = meanTraj(X,hittrials(slowTrials),6)'; %trajectory variable and predefined conditional trial indexes
rs = meanTraj(X,hittrials(fastTrials),6)'; %trajectory variable and predefined conditional trial indexes
%%% calculate initial condition differences in rf and rs
% rfrsinit = abs(rf(1,:)-rs(1,:));
% rs = rs-rfrrsinitl;
[neuralsimRT,neuraldiffRT,rprimef,rprimes] = neuralTrajDiff(rf,rs);
% reference to all hit trial initial state
[neuralTrajfh,rprimefinit,rprimehinit] = neuralTrajDiff(rf,rh,'initial');
[neuralTrajsh,rprimesinit,rprimehinit] = neuralTrajDiff(rs,rh,'initial');

%%% Now do analysis for MI hit vs FA
X = neuralTrajMIHitFA;
MIHittrials = 1:length(Behaviour.MIHitTrace);
MIFAtrials = length(Behaviour.MIHitTrace)+1:size(X,3);
rmih = meanTraj(X,MIHittrials,6)'; %trajectory variable and predefined conditional trial indexes
rmif = meanTraj(X,MIFAtrials,6)'; %trajectory variable and predefined conditional trial indexes
rmihfinit = abs(rmih(1,:)-rmif(1,:));
rmif = rmif-rmihfinit;
t = linspace(-Behaviour.parameters.windowBeforePull,Behaviour.parameters.windowAfterPull,size(rmih,1));
% Find important indices in array
stimStart = interp1(t,1:length(t),0,'nearest'); % Find zero of data which reflects some task condition we made
rawreactionTime = arrayfun(@(x) x.reactionTime,Behaviour.MIHitTrace);
mreactionTime = interp1(t,1:length(t),mean(rawreactionTime),'nearest');
reactionTime = interp1(t,1:length(t),rawreactionTime,'nearest');
reactionTime(isnan(reactionTime)) = length(t);
rawrewardTime = cell2mat(arrayfun(@(x) x.rewardIndex - x.LFPIndex(1501),Behaviour.MIHitTrace,'UniformOutput',false))/1000;
mrewardTime = interp1(t,1:length(t),mean(rawrewardTime),'nearest');
mrewardTime(isnan(reactionTime)) = length(t);


% Calculate differences in trajectories r'c(t)/||r'c(t)||
[neuralSimMI,neuralDiffMI,rprimemih,rprimemif] = neuralTrajDiff(rmih,rmif);

% reference to all hit trial initial state
[neuralTrajfh,rprimefinit,rprimehinit] = neuralTrajDiff(rf,rh,'initial');
[neuralTrajsh,rprimesinit,rprimehinit] = neuralTrajDiff(rs,rh,'initial');


showplot = 1;
if showplot
    %%% Hit vs Miss
    drawArrow = @(x,y,color) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'Color',color,'LineWidth',2,'MaxHeadSize',0.5);    
    figure,hold on
    clf
    % Updating the line
    x = rh(:,1);y = rh(:,2);z = rh(:,3);
    plot3(x,y,z,'-','color',[252/255 186/255 3/255],'lineWidth',1);hold on
    scatter3(x(stimStart,:),y(stimStart,:),z(stimStart,:),15,'g','filled')
    scatter3(x(mreactionTime),y(mreactionTime),z(mreactionTime),15,'b','filled')
    scatter3(x(mrewardTime),y(mrewardTime),z(mrewardTime),15,'r','filled')
    x = rm(:,1);y = rm(:,2);z = rm(:,3);
    plot3(x,y,z,'-','color',[3/255 190/255 252/255],'lineWidth',1);
    hold on,axis tight
    scatter3(x(1,:),y(1,:),z(1,:),15,'r','filled')
    scatter3(x(stimStart,:),y(stimStart,:),z(stimStart,:),15,'g','filled')
    %%% Hit vs FA
    figure
    clf
    % Updating the line
    x = rmih(:,1);y = rmih(:,2);z = rmih(:,3);
    plot3(x,y,z,'-','color',[252/255 186/255 3/255],'lineWidth',1);
    hold on,axis tight
    scatter3(x(stimStart,:),y(stimStart,:),z(stimStart,:),15,'g','filled')
    scatter3(x(mreactionTime),y(mreactionTime),z(mreactionTime),15,'b','filled')
    scatter3(x(mrewardTime),y(mrewardTime),z(mrewardTime),15,'r','filled')
    x = rmif(:,1);y = rmif(:,2);z = rmif(:,3);
    plot3(x,y,z,'-','color',[3/255 190/255 252/255],'lineWidth',1);
    hold on,axis tight
    scatter3(x(1,:),y(1,:),z(1,:),15,'r','filled')
    scatter3(x(stimStart,:),y(stimStart,:),z(stimStart,:),15,'g','filled')
    %     legend('Hit','','','','False Alarms')
    
    %%% Plot 2d PCA space with wave properties
    % Hit
    figure,subplot(131),plot(t,squeeze(neuralTrajHitMiss(1,:,hittrials)),'Color',[0 0 0 .25]);set(gca,'TickDir','out','fontsize',16),box off
    hold on, subplot(131),for n = 1:length(reactionTime),plot(t(reactionTime(n)),squeeze(neuralTrajHitMiss(1,reactionTime(n),n)),'b.'),end
    subplot(131),plot(t,mean(squeeze(neuralTrajHitMiss(1,:,hittrials)),2),'r','LineWidth',2)
    subplot(132),plot(t,squeeze(neuralTrajHitMiss(2,:,hittrials)),'Color',[0 0 0 .25]);set(gca,'TickDir','out','fontsize',16),box off
    hold on, subplot(132),for n = 1:length(reactionTime),plot(t(reactionTime(n)),squeeze(neuralTrajHitMiss(2,reactionTime(n),n)),'b.'),end
    subplot(132),plot(t,mean(squeeze(neuralTrajHitMiss(2,:,hittrials)),2),'r','LineWidth',2)
    subplot(133),plot(t,squeeze(neuralTrajHitMiss(3,:,hittrials)),'Color',[0 0 0 .25]);set(gca,'TickDir','out','fontsize',16),box off
    hold on, subplot(133),for n = 1:length(reactionTime),plot(t(reactionTime(n)),squeeze(neuralTrajHitMiss(3,reactionTime(n),n)),'b.'),end
    subplot(133),plot(t,mean(squeeze(neuralTrajHitMiss(3,:,hittrials)),2),'r','LineWidth',2)
    figure,subplot(131),plot(t,smoothdata(mean(waveDensityhit))/20*50,'k','LineWidth',2),set(gca,'TickDir','out','fontsize',16),box off
    subplot(131),hold on,plot(t,mean(waveDensityhit)/20*50,'b.')
    subplot(132),plot(t,smoothdata(mean(wavePGDhit)),'k','LineWidth',2),set(gca,'TickDir','out','fontsize',16),box off
    subplot(132),hold on,plot(t,mean(wavePGDhit),'b.')
    subplot(133),plot(t,smoothdata(mean(waveSpeedhit)),'k','LineWidth',2),set(gca,'TickDir','out','fontsize',16),box off
    subplot(133),hold on,plot(t,mean(waveSpeedhit),'b.')
    
    % Miss
    figure,subplot(131),plot(t,squeeze(neuralTrajHitMiss(1,:,misstrials)),'Color',[0 0 0 .25]);set(gca,'TickDir','out','fontsize',16),box off
    hold on,subplot(131),plot(t,mean(squeeze(neuralTrajHitMiss(1,:,misstrials)),2),'r','LineWidth',2)
    subplot(132),plot(t,squeeze(neuralTrajHitMiss(2,:,misstrials)),'Color',[0 0 0 .25]);set(gca,'TickDir','out','fontsize',16),box off
    hold on,subplot(132),plot(t,mean(squeeze(neuralTrajHitMiss(2,:,misstrials)),2),'r','LineWidth',2)
    subplot(133),plot(t,squeeze(neuralTrajHitMiss(3,:,misstrials)),'Color',[0 0 0 .25]);set(gca,'TickDir','out','fontsize',16),box off
    hold on,subplot(133),plot(t,mean(squeeze(neuralTrajHitMiss(3,:,misstrials)),2),'r','LineWidth',2)
    figure,subplot(131),plot(t,smoothdata(mean(waveDensitymiss))/20*50,'k','LineWidth',2),set(gca,'TickDir','out','fontsize',16),box off
    subplot(131),hold on,plot(t,mean(waveDensitymiss)/20*50,'b.')
    subplot(132),plot(t,smoothdata(mean(wavePGDmiss)),'k','LineWidth',2),set(gca,'TickDir','out','fontsize',16),box off
    subplot(132),hold on,plot(t,mean(wavePGDmiss),'b.')
    subplot(133),plot(t,smoothdata(mean(waveSpeedmiss)),'k','LineWidth',2),set(gca,'TickDir','out','fontsize',16),box off
    subplot(133),hold on,plot(t,mean(waveSpeedmiss),'b.')
end

%% Estimate differene in neural conditions
rprime1 = smoothdata(diff(squeeze(mean(neuralTrajHit(1,:,:),3)))-diff(squeeze(mean(neuralTrajMiss(1,:,:),3))),'gaussian',5);
rprime2 = smoothdata(diff(squeeze(mean(neuralTrajHit(2,:,:),3)))-diff(squeeze(mean(neuralTrajMiss(2,:,:),3))),'gaussian',5);
rprime3 = smoothdata(diff(squeeze(mean(neuralTrajHit(3,:,:),3)))-diff(squeeze(mean(neuralTrajMiss(3,:,:),3))),'gaussian',5);

figure,plot(-74*20:20:20*75,[0,rprime1],'LineWidth',2),hold on
plot(-74*20:20:20*75,[0,rprime2],'LineWidth',2)
plot(-74*20:20:20*75,[0,rprime3],'LineWidth',2),box off,set(gca,'FontSize',16),set(gca,'TickDir','out'),ylabel("r'"),


rprime1 = smoothdata(diff(squeeze(mean(neuralTrajMIHit(1,:,:),3)))-diff(squeeze(mean(neuralTrajMIFA(1,:,:),3))),'gaussian',5);
rprime2 = smoothdata(diff(squeeze(mean(neuralTrajMIHit(2,:,:),3)))-diff(squeeze(mean(neuralTrajMIFA(2,:,:),3))),'gaussian',5);
rprime3 = smoothdata(diff(squeeze(mean(neuralTrajMIHit(3,:,:),3)))-diff(squeeze(mean(neuralTrajMIFA(3,:,:),3))),'gaussian',5);

figure,plot(-74*20:20:20*75,[0,rprime1],'LineWidth',2),hold on
plot(-74*20:20:20*75,[0,rprime2],'LineWidth',2)
plot(-74*20:20:20*75,[0,rprime3],'LineWidth',2),box off,set(gca,'FontSize',16),set(gca,'TickDir','out'),ylabel("r'")
%% Overlaying traveling wave dynamics across neural trajectories
% TODO: Plotting the data like this makes the rendering all messed up; need
% to adapt from Lyles GP phase code for plotting....
figure
x = rh(:,1);
y = rh(:,2);
p = plot(x,y,'r', 'LineWidth',2);
% modified jet-colormap
% cd = [uint8(jet(150)*255) uint8(ones(150,1))].';
drawnow
set(p.Edge, 'ColorBinding','interpolated', 'ColorData',uint8(cd1))
%% index PGD
cd = [uint8(jet(150)*255) uint8(ones(150,1))].';
n = 150;
col = mean(waveDensityhit)';
col_map = (col - min(col))/(max(col)-min(col)) * (n-1) + 1;
for n = 1:150
    cd1(:,n) = cd(:,floor(col_map(n)));
end

end




%% Local functions

function r = meanTraj(X,trials,components)
r = squeeze(mean(X(1:components,:,trials),3));
end

function [neuralTrajSim,neuralTrajdiff,rprimehnorm,rprimemnorm] = neuralTrajDiff(r1,r2,varargin)
if strcmp(varargin,'initial'),initCalc = 1;else, initCalc = 0;end
if initCalc %r'(t)/||r'(t)|| * r(0)-ri(t)/||r(0)-ri(t)||
    rprimeh = diff(r1);
    rprimehnorm = rprimeh/norm(rprimeh);
    rprimem = r2(1,:)-r1; % control condition for reference
    rprimemnorm = rprimem/norm(rprimem);
    neuralTrajSim = dot(rprimehnorm',rprimemnorm(1:end-1,:)');
    neuralTrajdiff = rprimeh-rprimem;
else
    rprimeh = diff(r1);
    rprimehnorm = rprimeh/norm(rprimeh);
    rprimem = diff(r2);
    rprimemnorm = rprimem/norm(rprimem);
    neuralTrajSim = dot(rprimehnorm',rprimemnorm');
    neuralTrajdiff = rprimeh-rprimem;
end
end