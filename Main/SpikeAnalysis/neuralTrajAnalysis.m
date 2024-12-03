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

% X = arrayfun(@(x) vertcat(x.xorth),Spikes.GPFA.seqTrainBaselineOpto,'UniformOutput',false);
% X = horzcat(X{:});
% neuralTrajBaselineOpto = reshape(X,size(X,1),Spikes.GPFA.seqTrainBaselineOpto(1).T,[]);

%% Calculate average trajectories and divergence based on trial difference
% local function call for meaning based on combined PCA of trial conditions
X = neuralTrajHitMiss;
hittrials = 1:length(Behaviour.cueHitTrace);
misstrials = length(Behaviour.cueHitTrace)+1:size(X,3);
rh = meanTraj(X,hittrials,6)'; %trajectory variable and predefined conditional trial indexes
rm = meanTraj(X,misstrials,6)'; %trajectory variable and predefined conditional trial indexes
%%
[r,s] = meanTraj(X,hittrials,6);
%%
% rhminit = abs(rh(1,:)-rm(1,:));
% rm = rm-rhminit;
t = linspace(-Behaviour.parameters.windowBeforeCue,Behaviour.parameters.windowAfterCue,size(rh,1));
% Find important indices in array
stimStart = interp1(t,1:length(t),0,'nearest'); % Find zero of data which reflects some task condition we made
rawreactionTime = arrayfun(@(x) x.reactionTime,Behaviour.cueHitTrace);
mreactionTime = interp1(t,1:length(t),mean(rawreactionTime),'nearest');
reactionTime = interp1(t,1:length(t),rawreactionTime,'nearest');
reactionTime(isnan(reactionTime)) = length(t);
if ~isempty(Waves)
    rawWaveDensityhit = arrayfun(@(x) vertcat(x.wavePresent),Waves.wavesHit,'UniformOutput',false);rawWaveDensityhit= vertcat(rawWaveDensityhit{:});
    rawWavePGDhit = arrayfun(@(x) vertcat(x.PGD), Waves.wavesHit,'UniformOutput',false);rawWavePGDhit = vertcat(rawWavePGDhit{:});
    rawWaveSpeedhit = arrayfun(@(x) vertcat(x.s), Waves.wavesHit,'UniformOutput',false);rawWaveSpeedhit = vertcat(rawWaveSpeedhit{:});
    
    rawWaveDensitymiss = arrayfun(@(x) vertcat(x.wavePresent),Waves.wavesMiss,'UniformOutput',false);rawWaveDensitymiss= vertcat(rawWaveDensitymiss{:});
    rawWavePGDmiss = arrayfun(@(x) vertcat(x.PGD), Waves.wavesMiss,'UniformOutput',false);rawWavePGDmiss = vertcat(rawWavePGDmiss{:});
    rawWaveSpeedmiss = arrayfun(@(x) vertcat(x.s), Waves.wavesMiss,'UniformOutput',false);rawWaveSpeedmiss = vertcat(rawWaveSpeedmiss{:});
    
    rawWaveDensityFA = arrayfun(@(x) vertcat(x.wavePresent),Waves.wavesMIFA,'UniformOutput',false);rawWaveDensityFA= vertcat(rawWaveDensityFA{:});
    rawWavePGDFA = arrayfun(@(x) vertcat(x.PGD), Waves.wavesMIFA,'UniformOutput',false);rawWavePGDFA = vertcat(rawWavePGDFA{:});
    rawWaveSpeedFA = arrayfun(@(x) vertcat(x.s), Waves.wavesMIFA,'UniformOutput',false);rawWaveSpeedFA = vertcat(rawWaveSpeedFA{:});
    
    rawWaveDensityMIHit = arrayfun(@(x) vertcat(x.wavePresent),Waves.wavesMIHit,'UniformOutput',false);rawWaveDensityMIHit= vertcat(rawWaveDensityMIHit{:});
    rawWavePGDMIHit= arrayfun(@(x) vertcat(x.PGD), Waves.wavesMIHit,'UniformOutput',false);rawWavePGDMIHit = vertcat(rawWavePGDMIHit{:});
    rawWaveSpeedMIHit = arrayfun(@(x) vertcat(x.s), Waves.wavesMIHit,'UniformOutput',false);rawWaveSpeedMIHit = vertcat(rawWaveSpeedMIHit{:});
    
    
    waveDensityhit = [];
    wavePGDhit = [];
    waveSpeedhit = [];
    
    waveDensitymiss = [];
    wavePGDmiss = [];
    waveSpeedmiss = [];
    
    waveDensityFA = [];
    wavePGDFA = [];
    waveSpeedFA = [];
    
    waveDensityMIHit = [];
    wavePGDMIHit = [];
    waveSpeedMIHit = [];
    
    win = ceil(1:20:size(rawWaveDensityhit,2));
    for n = 1:length(win)-1
        waveDensityhit = horzcat(waveDensityhit,sum(rawWaveDensityhit(:,win(n):win(n+1)),2)); %convert to wave/sec
        wavePGDhit = horzcat(wavePGDhit,mean(rawWavePGDhit(:,win(n):win(n+1)),2));
        waveSpeedhit = horzcat(waveSpeedhit,mean(rawWaveSpeedhit(:,win(n):win(n+1)),2));
        
        waveDensitymiss = horzcat(waveDensitymiss,sum(rawWaveDensitymiss(:,win(n):win(n+1)),2)); %convert to wave/sec
        wavePGDmiss = horzcat(wavePGDmiss,mean(rawWavePGDmiss(:,win(n):win(n+1)),2));
        waveSpeedmiss = horzcat(waveSpeedmiss,mean(rawWaveSpeedmiss(:,win(n):win(n+1)),2));
        
        waveDensityFA = horzcat(waveDensityFA,sum(rawWaveDensityFA(:,win(n):win(n+1)),2)); %convert to wave/sec
        wavePGDFA = horzcat(wavePGDFA,mean(rawWavePGDFA(:,win(n):win(n+1)),2));
        waveSpeedFA = horzcat(waveSpeedFA,mean(rawWaveSpeedFA(:,win(n):win(n+1)),2));
        
        waveDensityMIHit = horzcat(waveDensityMIHit,sum(rawWaveDensityMIHit(:,win(n):win(n+1)),2)); %convert to wave/sec
        wavePGDMIHit = horzcat(wavePGDMIHit,mean(rawWavePGDMIHit(:,win(n):win(n+1)),2));
        waveSpeedMIHit = horzcat(waveSpeedMIHit,mean(rawWaveSpeedMIHit(:,win(n):win(n+1)),2));
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

% rb = meanTraj(X,hittrials(Behaviour.hitTemp>-10),6)';
% rc = meanTraj(X,hittrials(Behaviour.hitTemp<=-10),6)';

%%% calculate initial condition differences in rf and rs
% rfrsinit = abs(rf(1,:)-rs(1,:));
% rs = rs-rfrrsinitl;
[neuralsimRT,neuraldiffRT,rprimef,rprimes] = neuralTrajDiff(rf,rs);
%%
% [neuralsimRT,neuraldiffRT,rprimef,rprimes] = neuralTrajDiff(rb,rc);
%%

% reference to all hit trial initial state
% [neuralTrajfh,rprimefinit,rprimehinit] = neuralTrajDiff(rf,rh,'initial');
% [neuralTrajsh,rprimesinit,rprimehinit] = neuralTrajDiff(rs,rh,'initial');

%%% Now do analysis for MI hit vs FA
X = neuralTrajMIHitFA;
MIHittrials = 1:length(Behaviour.MIHitTrace);
MIFAtrials = length(Behaviour.MIHitTrace)+1:size(X,3);
rmih = meanTraj(X,MIHittrials,6)'; %trajectory variable and predefined conditional trial indexes
rmif = meanTraj(X,MIFAtrials,6)'; %trajectory variable and predefined conditional trial indexes
%%
[r,s] = meanTraj(X,MIFAtrials,6);
%%
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

x = rmif(:,1);
y = rmif(:,2);
cd = [uint8((jet(150))*255) uint8(ones(150,1))].';
n = 150;
t = mean(wavePGDhit);
col = [smoothdata(diff(t),'movmean',10)'];
%col = smoothdata(t,'movmean',10);
% Interp to make the line smoother
xin = interp1(1:150,x,1:0.1:150);
yin = interp1(1:150,y,1:0.1:150);
col = interp1(1:150,col,1:0.1:150);
col_map = (col - min(col))/(max(col)-min(col)) * (n-1) + 1;
for n = 1:length(col)
    cd1(:,n) = cd(:,floor(col_map(n)));
end
figure,
for n = 2:length(col)
plot(xin(n-1:n),yin(n-1:n),'color',double(cd1(1:3,n))/255, 'LineWidth',2);hold on %cline( time, xf, [], angle(xgp) );
end
box off, axis off
% modified jet-colormap
% cd = [uint8(jet(150)*255) uint8(ones(150,1))].';
%% Use Cline so we can make the colorbar
load myMap
figure,
h4 = cline( xin, yin, [], col);
colormap((jet))
set( h4, 'linestyle', '-', 'linewidth', 2  );axis off

%% Statistics of Neural Traj and Waves
X = neuralTrajHitMiss;
PQ = [];
idx = discretize(reactionTime,20);

for n = 1:length(hittrials)
dat = meanTraj(X,n,6); %grab and calculate distance per reaction time
x = dat(1,1:reactionTime(n)); y = dat(2,1:reactionTime(n)); z = dat(3,1:reactionTime(n));
PQ(n) = sqrt((x(end)-x(1))^2+(y(end)-y(1))^2+(z(end)-z(1))^2); %Calculate Speed
dat1 = vertcat(zeros(1,6),rprimeh)'; % Grab norm vector alignment data
x = dat1(1,1:reactionTime(n)); y = dat1(2,1:reactionTime(n)); z = dat1(3,1:reactionTime(n));
TrajAngle(n) = mean(mean([x;y;z])); % Calculates how well aligned trajectories are for all trials
if ~isempty(Waves)
    dPGD(n) = mean(diff(wavePGDhit(n,stimStart:reactionTime(n))));
end
end
%%
X = neuralTrajHitMiss;
PQ = [];
idx = discretize(reactionTime,20);

for n = 1:length(hittrials)
dat = meanTraj(X,n,6); %grab and calculate distance per reaction time
x = dat(1,:); y = dat(2,:); z = dat(3,:);
PQ(n) = sqrt((x(end)-x(1))^2+(y(end)-y(1))^2+(z(end)-z(1))^2); %Calculate Distance
dat1 = vertcat(zeros(1,6),rprimeh)'; % Grab norm vector alignment data
x = dat1(1,1:reactionTime(n)); y = dat1(2,1:reactionTime(n)); z = dat1(3,1:reactionTime(n));
TrajAngle(n) = mean(mean([x;y;z])); % Calculates how well aligned trajectories are for all trials
if ~isempty(Waves)
    dPGD(n) = mean(diff(wavePGDhit(n,stimStart:reactionTime(n))));
end
end
%%
dat1 = PQ(Behaviour.hitTemp>-10);
dat2 = PQ(Behaviour.hitTemp<-10);
temp = nan(max([length(dat1) length(dat2)]),2);
temp(1:length(dat1),1) = dat1;
temp(1:length(dat2),2) = dat2;

figure,customBoxplot(temp)
%%
PQ = (PQ./rawreactionTime)';

figure,scatter(TrajAngle,rawreactionTime,'k','filled')
xlim([0 .02])
ylim([0 2])
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
xlabel('Neural trajectory alignment'),ylabel('Reaction time (s)')

figure,scatter(PQ,rawreactionTime,'k','filled')
xlim([0 15])
ylim([0 2])
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
xlabel('Neural trajectory speed'),ylabel('Reaction time (s)')

mdl = fitlm(PQ,rawreactionTime)
figure,plot(mdl)
xlim([0 15])
ylim([0 2])
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
xlabel('Neural trajectory speed'),ylabel('Reaction time (s)')
legend(num2str(mdl.Rsquared.Ordinary),num2str(mdl.Coefficients.pValue(2)))
title('')
%%
figure,scatter(dPGD,TrajAngle,'k','filled')
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
xlabel('Neural trajectory speed'),ylabel('Phase gradient')
xlim([0 10])

figure,scatter(PQ,dPGD,'k','filled')
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
xlabel('Neural trajectory speed'),ylabel('Phase gradient')
xlim([0 10])

mdl = fitlm(PQ,dPGD)
figure,plot(mdl)
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
xlabel('Neural trajectory speed'),ylabel('Phase gradient')
legend(num2str(mdl.Rsquared.Ordinary),num2str(mdl.Coefficients.pValue(2)))
xlim([0 10])
title('')

figure,scatter(dPGD,rawreactionTime,'k','filled')
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
xlabel('Phase gradient'),ylabel('Reaction time (s)')

mdl = fitlm(dPGD,rawreactionTime)
figure,plot(mdl)
set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
xlabel('Phase gradient'),ylabel('Reaction time (s)')
legend(num2str(mdl.Rsquared.Ordinary),num2str(mdl.Coefficients.pValue(2)))
%%
twoD = 1;
showplot = 1;
colors = [0 0.4470 0.7410;0.75 0.75 0.75;190/255 30/255 45/255];
if showplot
    if twoD
        %%% Hit vs Miss
        drawArrow = @(x,y,color) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'Color',color,'LineWidth',2,'MaxHeadSize',0.5);
        figure,hold on
        clf
        % Updating the line
        x = rh(:,1);y = rh(:,2);z = rh(:,3);
        plot(x,y,'-','color',colors(1,:),'lineWidth',2);hold on
        scatter(x(stimStart,:),y(stimStart,:),25,'k','filled')
        %scatter3(x(mreactionTime),y(mreactionTime),z(mreactionTime),15,'b','filled')
        %scatter3(x(mrewardTime),y(mrewardTime),z(mrewardTime),15,'r','filled')
        %scatter3(x(1:5:150,:),y(1:5:150,:),z(1:5:150,:),'k')
        x = rm(:,1);y = rm(:,2);z = rm(:,3);
        plot(x,y,'-','color',colors(2,:),'lineWidth',2);
        hold on,axis tight
        scatter3(x(stimStart,:),y(stimStart,:),25,'k','filled')
        set(gca,'tickdir','out')
        %scatter3(x(1:5:150,:),y(1:5:150,:),z(1:5:150,:),'k')
        %scatter3(x(mreactionTime,:),y(mreactionTime,:),z(mreactionTime,:),15,'g','filled')
        %%% Hit vs FA
        figure
        clf
        % Updating the line
        x = rmih(:,2);y = rmih(:,3);z = rmih(:,3);
        plot(x,y,'-','color',colors(1,:),'lineWidth',2);
        hold on,axis tight
        scatter(x(stimStart,:),y(stimStart,:),25,'k','filled')
        %scatter3(x(mreactionTime),y(mreactionTime),z(mreactionTime),15,'b','filled')
        %scatter3(x(mrewardTime),y(mrewardTime),z(mrewardTime),15,'r','filled')
        %scatter3(x(1:5:150,:),y(1:5:150,:),z(1:5:150,:),'k')
        x = rmif(:,1);y = rmif(:,2);z = rmif(:,3);
        plot(x,y,'-','color',colors(3,:),'lineWidth',2);
        hold on,axis tight
        scatter(x(1,:),y(1,:),25,'k','filled')
        scatter(x(stimStart,:),y(stimStart,:),25,'k','filled')
        
    else
        %%% Hit vs Miss
        drawArrow = @(x,y,color) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'Color',color,'LineWidth',2,'MaxHeadSize',0.5);
        figure,hold on
        clf
        % Updating the line
        x = rh(:,1);y = rh(:,2);z = rh(:,3);
        plot3(x,y,z,'-','color',colors(1,:),'lineWidth',2);hold on
        scatter3(x(stimStart,:),y(stimStart,:),z(stimStart,:),15,'k','filled')
        %scatter3(x(mreactionTime),y(mreactionTime),z(mreactionTime),15,'k','filled')
        %scatter3(x(mrewardTime),y(mrewardTime),z(mrewardTime),15,'k','filled')
        %scatter3(x(1:5:150,:),y(1:5:150,:),z(1:5:150,:),'k')
        x = rm(:,1);y = rm(:,2);z = rm(:,3);
        plot3(x,y,z,'-','color',colors(2,:),'lineWidth',2);
        hold on,axis tight
        scatter3(x(stimStart,:),y(1,:),z(1,:),15,'k','filled')
        %scatter3(x(mreactionTime,:),y(mreactionTime,:),z(mreactionTime,:),15,'k','filled')
        view([0,0])
         set(gca,'tickdir','out'),axis square
        %%% Hit vs FA
        figure
        clf
        % Updating the line
        x = rmih(:,1);y = rmih(:,2);z = rmih(:,3);
        plot3(x,y,z,'-','color',colors(1,:),'lineWidth',2);
        hold on,axis tight
        scatter3(x(stimStart,:),y(stimStart,:),z(stimStart,:),15,'k','filled')
        %scatter3(x(mreactionTime),y(mreactionTime),z(mreactionTime),15,'b','filled')
        %scatter3(x(mrewardTime),y(mrewardTime),z(mrewardTime),15,'r','filled')
        %scatter3(x(1:5:150,:),y(1:5:150,:),z(1:5:150,:),'k')
        x = rmif(:,1);y = rmif(:,2);z = rmif(:,3);
        plot3(x,y,z,'-','color',colors(3,:),'lineWidth',2);
        hold on,axis tight
        scatter3(x(1,:),y(1,:),z(1,:),15,'r','filled')
        scatter3(x(stimStart,:),y(stimStart,:),z(stimStart,:),15,'g','filled')
        %     legend('Hit','','','','False Alarms')
        view([0,0])
         set(gca,'tickdir','out'),axis square
    end
    %% Plot 2d PCA space with wave properties
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
end




%% Local functions

function [r,s] = meanTraj(X,trials,components)
matrix = X(1:components,:,trials);
r = squeeze(mean(matrix,3));
% Subtract the mean from each element
centered_matrix = matrix - r;

% Calculate the Euclidean distance for each row
s = squeeze(sqrt(sum(centered_matrix.^2, 2)));
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