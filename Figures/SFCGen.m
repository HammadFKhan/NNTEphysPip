%% GP-SFC
if ~exist('SFC','var')
    SFC = struct();
    files = dir(fullfile('D:\M1_GSP\','*.mat'));
    for fileNum = 1:length(files)
        fName = fullfile(files(fileNum).folder,files(fileNum).name);
        disp(['Loading ' fName '...'])
        load(fName)
        % Extract GP-SFC values and make combined structure for analysis
        SFC(fileNum).name = files(fileNum).name;
        SFC(fileNum).nSpikes = Spikes.nSpikes;
        SFC(fileNum).SpikeField = Spikes.SpikeField;
    end
end

%% Calculate SPI across electrodes (maybe we show as a function of max SPI
%%% per neurons?
figure,imagesc(SFC(sess).SpikeField.postcuehit.SPI)
%% Stats of SPI and angle
% SPI and angle of all neurons relative to their closest electrodes
Stats.hit.preStimSPI = [];Stats.hit.preStimAng = [];Stats.hit.postStimSPI = [];Stats.hit.postStimAng = [];
Stats.miss.preStimSPI = [];Stats.miss.preStimAng = [];Stats.miss.postStimSPI = [];Stats.miss.postStimAng = [];
Stats.FA.preStimSPI = [];Stats.FA.preStimAng = [];Stats.FA.postStimSPI = [];Stats.FA.postStimAng = [];

for sess = 1:9
    pos = SFC(sess).SpikeField.precuehit.SpikeElectrode;
    if max(pos)>size(SFC(sess).SpikeField.precuehit.SPI,1)
        pos = ceil(pos/2);
    end
    for n = 1:SFC(sess).nSpikes
        Stats.hit.preStimSPI = [Stats.hit.preStimSPI;SFC(sess).SpikeField.precuehit.SPI(pos(n),n)];
        Stats.hit.preStimAng = [Stats.hit.preStimAng;SFC(sess).SpikeField.precuehit.Angle(pos(n),n)];
        Stats.hit.postStimSPI = [Stats.hit.postStimSPI;max(SFC(sess).SpikeField.postcuehit.SPI(:,n))];
        Stats.hit.postStimAng = [Stats.hit.postStimAng;SFC(sess).SpikeField.postcuehit.Angle(pos(n),n)];
        
        Stats.miss.preStimSPI = [Stats.miss.preStimSPI;SFC(sess).SpikeField.precuemiss.SPI(pos(n),n)];
        Stats.miss.preStimAng = [Stats.miss.preStimAng;SFC(sess).SpikeField.precuemiss.Angle(pos(n),n)];
        Stats.miss.postStimSPI = [Stats.miss.postStimSPI;max(SFC(sess).SpikeField.postcuemiss.SPI(:,n))];
        Stats.miss.postStimAng = [Stats.miss.postStimAng;max(SFC(sess).SpikeField.postcuemiss.Angle(:,n))];
        
        Stats.FA.preStimSPI = [Stats.FA.preStimSPI;SFC(sess).SpikeField.precueMIFA.SPI(pos(n),n)];
        Stats.FA.preStimAng = [Stats.FA.preStimAng;SFC(sess).SpikeField.precueMIFA.Angle(pos(n),n)];
        Stats.FA.postStimSPI = [Stats.FA.postStimSPI;SFC(sess).SpikeField.postcueMIFA.SPI(pos(n),n)];
        Stats.FA.postStimAng = [Stats.FA.postStimAng;SFC(sess).SpikeField.postcueMIFA.Angle(pos(n),n)];
    end
    %  figure(1)
%     for n = 1:SFC(sess).nSpikes
%         scatter(dat(pos(n),n),pos(n),sz(pos(n),n)*150,'b','filled'),hold on
%         set(gca, 'YDir','reverse')
%     end
end
%%
colors = [0 0.4470 0.7410;0.75 0.75 0.75;190/255 30/255 45/255];

temp = [Stats.hit.postStimSPI Stats.miss.postStimSPI Stats.FA.postStimSPI];
figure(1)
clf
violinplot(temp,[],'ShowData',true,'ShowWhiskers',false,'ShowBox',false,'MarkerSize',5,'ViolinColor',colors);
set(gca,'tickdir','out','fontsize',16),box off,axis square
ylabel('Spike Field Coherence'),ylim([0 1])
[p,~,stats] = anova1(temp);
%%
temp = [Stats.hit.postStimSPI-Stats.hit.preStimSPI  Stats.miss.postStimSPI-Stats.miss.preStimSPI Stats.FA.postStimSPI-Stats.FA.preStimSPI];
figure(2)
clf
violinplot(temp,[],'ShowData',true,'ShowWhiskers',false,'ShowBox',false,'MarkerSize',5,'ViolinColor',colors);
set(gca,'tickdir','out','fontsize',16),box off
ylabel('Post-Pre Stimulus SFC')
ylim([-1 1])
axis square
p = anova1(temp);
%%
figure,histogram(preCueAng,-pi:pi/6:pi,'normalization','probability','edgecolor','none'),hold on
histogram(postCueAng,-pi:pi/6:pi,'normalization','probability','edgecolor','none')

%% Here we see how does the spiking phase and amplitude change pre- post- cue

spiDiff = [];
angDiff = [];
for sess = 1:9
    spiPre = SFC(sess).SpikeField.precuehit.SPI;
    angPre = SFC(sess).SpikeField.precuehit.Angle;
    spiPost = SFC(sess).SpikeField.postcuehit.SPI;
    angPost = SFC(sess).SpikeField.postcuehit.Angle;
    pos = SFC(sess).SpikeField.precuehit.SpikeElectrode;
    if max(pos)>size(dat,1)
        pos = ceil(pos/2);
    end
    for n = 1:SFC(sess).nSpikes 
    spiDiff = [spiDiff, spiPost(pos(n),n)-spiPre(pos(n),n)];
    angDiff = [angDiff, angPost(pos(n),n)-angPre(pos(n),n)];
    end
end

figure,histogram(spiDiff,10,'normalization','probability','edgecolor','none')
figure,histogram(angDiff,-pi:pi/6:pi,'normalization','probability','edgecolor','none')
%% Plot example GP-SFC
if 1==2
    load('D:\M1LFPCued\23025_Rbp4_Day3LFP.mat') % Selectively load this
    load('D:\M1_GSP\23025Rbp4_Day3_SpikesSF.mat')
end
figure(1),clf
imagesc(interp2(Spikes.SpikeField.postcuehit.PhaseCorr,1)) ,axis square,xlabel('Electrode'),ylabel('Electrode'),colorbar
set(gca,'fontsize',16),colormap(cool)
%%
f = figure(2);clf
plot(mean(squeeze(LFP.probe1.hitLFP{1}),1));
set(gca,'tickdir','out','fontsize',16),box off
f.Position = [681 559 560 220];
%%
trial = 52;
eld = 16;
time = 1:3001;
xw = squeeze(LFP.probe1.hitLFP{trial}(eld,:,:));
xgp = squeeze( LFP.probe1.hitxgp{trial}(eld,:,:)) ;
xf = squeeze(LFP.probe1.hitLFP{trial}(eld,:,:));

cd = [uint8((cool(3001))*255) uint8(ones(3001,1))].';
n = 3001;
col = [angle(xgp)];
%col = smoothdata(t,'movmean',10);
% Interp to make the line smoother
col_map = (col - min(col))/(max(col)-min(col)) * (n-1) + 1;
for n = 1:length(col)
    cd1(:,n) = cd(:,floor(col_map(n)));
end
f = figure;
for n = 2:length(col)
plot(time(n-1:n),xw(n-1:n),'color',double(cd1(1:3,n))/255, 'LineWidth',1);hold on %cline( time, xf, [], angle(xgp) );
end
box off, axis on
f.Position = [681 559 760 220];
xlim([1200 2400])
% inset
map = colorcet( 'C2' ); colormap( circshift( 'cool', [ 28, 0 ] ) )
ax2 = axes; set( ax2, 'position', [0.2116    0.6976    0.0884    0.2000] ); axis image
[x1,y1] = pol2cart( angle( exp(1i.*linspace(-pi,pi,100)) ), ones( 1, 100 ) );
h3 = cline( x1, y1, linspace(-pi,pi,100) ); axis off; set( h3, 'linewidth', 6 )

% text labels
t1 = text( 0, 0, 'GP' );
set( t1, 'fontname', 'arial', 'fontsize', 28, 'fontweight', 'bold', 'horizontalalignment', 'center' )
set( gcf, 'currentaxes', ax1 )


%% Make Spike population map
Spikes = makeSpikeGPFA(Spikes);
%%
f = figure(3);clf
Show_Spikes(Spikes.GPFA.hit.dat(trial).spikes);
set(gca,'tickdir','out','fontsize',16),box off
f.Position = [681 559 760 220];
xlim([1200 2400])
%% 
f = figure(4);clf
plot(smoothdata(IntanBehaviour.cueHitTrace(trial).trace)),axis off
f.Position = [681 559 760 220];
xlim([1200 2400])