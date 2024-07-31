%% GP-SFC
if ~exist('SFC','var')
    SFC = struct();
    files = dir(fullfile('D:\M1Cooling\SpikesGSP\','*.mat'));
    for fileNum = 1:length(files)
        fName = fullfile(files(fileNum).folder,files(fileNum).name);
        disp(['Loading ' fName '...'])
        load(fName)
        % Extract GP-SFC values and make combined structure for analysis
        SFC(fileNum).name = files(fileNum).name;
        SFC(fileNum).nSpikes = Spikes.nSpikes;
        SFC(fileNum).BaselineSpikeField = Spikes.Baseline.SpikeField;
        SFC(fileNum).CooledSpikeField = Spikes.Cooled.SpikeField;
        
    end
end

%% Calculate SPI across electrodes (maybe we show as a function of max SPI
%%% per neurons?
figure,imagesc(SFC(sess).CooledSpikeField.postcuehit.SPI)
%% Stats of SPI and angle
% SPI and angle of all neurons relative to their closest electrodes

Stats = sessionStats(SFC(1:2));
%% Plot violins of SFC coherence for baseline and cooled with anova statistics
colors = [0 0.4470 0.7410;0.75 0.75 0.75;190/255 30/255 45/255];
dat = Stats.baseline;
temp = [dat.hit.postStimSPI dat.miss.postStimSPI dat.FA.postStimSPI];
figure(1)
clf
violinplot(temp,[],'ShowData',true,'ShowWhiskers',false,'ShowBox',false,'MarkerSize',5,'ViolinColor',colors);
set(gca,'tickdir','out','fontsize',16),box off,axis square
ylabel('Spike Field Coherence'),ylim([0 1])
[p,~,stats] = anova1(temp);

colors = [0 0.4470 0.7410;0.75 0.75 0.75;190/255 30/255 45/255];
dat = Stats.cooled;
temp = [dat.hit.postStimSPI dat.miss.postStimSPI dat.FA.postStimSPI];
figure(10)
clf
violinplot(temp,[],'ShowData',true,'ShowWhiskers',false,'ShowBox',false,'MarkerSize',5,'ViolinColor',colors);
set(gca,'tickdir','out','fontsize',16),box off,axis square
ylabel('Spike Field Coherence'),ylim([0 1])
[p,~,stats] = anova1(temp);
%%
close all
dat = Stats.baseline;
temp = [dat.hit.postStimSPI-dat.hit.preStimSPI  dat.miss.postStimSPI-dat.miss.preStimSPI dat.FA.postStimSPI-dat.FA.preStimSPI];
figure(1)
clf
violinplot(temp,[],'ShowData',true,'ShowWhiskers',false,'ShowBox',false,'MarkerSize',5,'ViolinColor',colors);
set(gca,'tickdir','out','fontsize',16),box off
ylabel('Post-Pre Stimulus SFC')
ylim([-1 1])
axis square
p = anova1(temp);

dat = Stats.cooled;
temp = [dat.hit.postStimSPI-dat.hit.preStimSPI  dat.miss.postStimSPI-dat.miss.preStimSPI dat.FA.postStimSPI-dat.FA.preStimSPI];
figure(10)
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
%% LOCAL FUNCTIONS

function Stats = sessionStats(SFC)

% Someone needs to revoke my MATLAB access for how bad this function is
% built



Stats.baseline.hit.preStimSPI = [];Stats.baseline.hit.preStimAng = [];Stats.baseline.hit.postStimSPI = [];Stats.baseline.hit.postStimAng = [];
Stats.baseline.miss.preStimSPI = [];Stats.baseline.miss.preStimAng = [];Stats.baseline.miss.postStimSPI = [];Stats.baseline.miss.postStimAng = [];
Stats.baseline.FA.preStimSPI = [];Stats.baseline.FA.preStimAng = [];Stats.baseline.FA.postStimSPI = [];Stats.baseline.FA.postStimAng = [];

Stats.cooled.hit.preStimSPI = [];Stats.cooled.hit.preStimAng = [];Stats.cooled.hit.postStimSPI = [];Stats.cooled.hit.postStimAng = [];
Stats.cooled.miss.preStimSPI = [];Stats.cooled.miss.preStimAng = [];Stats.cooled.miss.postStimSPI = [];Stats.cooled.miss.postStimAng = [];
Stats.cooled.FA.preStimSPI = [];Stats.cooled.FA.preStimAng = [];Stats.cooled.FA.postStimSPI = [];Stats.cooled.FA.postStimAng = [];

for sess = 1:length(SFC)
    spikeField = SFC(sess).BaselineSpikeField;
    pos = spikeField.precuehit.SpikeElectrode;
    if max(pos)>size(spikeField.precuehit.SPI,1)
        pos = ceil(pos/2);
    end
    for n = 1:SFC(sess).nSpikes
        Stats.baseline.hit.preStimSPI = [Stats.baseline.hit.preStimSPI;spikeField.precuehit.SPI(pos(n),n)];
        Stats.baseline.hit.preStimAng = [Stats.baseline.hit.preStimAng;spikeField.precuehit.Angle(pos(n),n)];
        Stats.baseline.hit.postStimSPI = [Stats.baseline.hit.postStimSPI;max(spikeField.postcuehit.SPI(:,n))];
        Stats.baseline.hit.postStimAng = [Stats.baseline.hit.postStimAng;spikeField.postcuehit.Angle(pos(n),n)];
        
        Stats.baseline.miss.preStimSPI = [Stats.baseline.miss.preStimSPI;spikeField.precuemiss.SPI(pos(n),n)];
        Stats.baseline.miss.preStimAng = [Stats.baseline.miss.preStimAng;spikeField.precuemiss.Angle(pos(n),n)];
        Stats.baseline.miss.postStimSPI = [Stats.baseline.miss.postStimSPI;max(spikeField.postcuemiss.SPI(:,n))];
        Stats.baseline.miss.postStimAng = [Stats.baseline.miss.postStimAng;spikeField.postcuemiss.Angle(pos(n),n)];
        
        Stats.baseline.FA.preStimSPI = [Stats.baseline.FA.preStimSPI;spikeField.precueMIFA.SPI(pos(n),n)];
        Stats.baseline.FA.preStimAng = [Stats.baseline.FA.preStimAng;spikeField.precueMIFA.Angle(pos(n),n)];
        Stats.baseline.FA.postStimSPI = [Stats.baseline.FA.postStimSPI;max(spikeField.postcueMIFA.SPI(:,n))];
        Stats.baseline.FA.postStimAng = [Stats.baseline.FA.postStimAng;spikeField.postcueMIFA.Angle(pos(n),n)];
    end
    
    spikeField = SFC(sess).CooledSpikeField;
    pos = spikeField.precuehit.SpikeElectrode;
    if max(pos)>size(spikeField.precuehit.SPI,1)
        pos = ceil(pos/2);
    end
    for n = 1:SFC(sess).nSpikes
        Stats.cooled.hit.preStimSPI = [Stats.cooled.hit.preStimSPI;spikeField.precuehit.SPI(pos(n),n)];
        Stats.cooled.hit.preStimAng = [Stats.cooled.hit.preStimAng;spikeField.precuehit.Angle(pos(n),n)];
        Stats.cooled.hit.postStimSPI = [Stats.cooled.hit.postStimSPI;max(spikeField.postcuehit.SPI(:,n))];
        Stats.cooled.hit.postStimAng = [Stats.cooled.hit.postStimAng;spikeField.postcuehit.Angle(pos(n),n)];
        
        Stats.cooled.miss.preStimSPI = [Stats.cooled.miss.preStimSPI;spikeField.precuemiss.SPI(pos(n),n)];
        Stats.cooled.miss.preStimAng = [Stats.cooled.miss.preStimAng;spikeField.precuemiss.Angle(pos(n),n)];
        Stats.cooled.miss.postStimSPI = [Stats.cooled.miss.postStimSPI;max(spikeField.postcuemiss.SPI(:,n))];
        Stats.cooled.miss.postStimAng = [Stats.cooled.miss.postStimAng;spikeField.postcuemiss.Angle(pos(n),n)];
        
        Stats.cooled.FA.preStimSPI = [Stats.cooled.FA.preStimSPI;spikeField.precueMIFA.SPI(pos(n),n)];
        Stats.cooled.FA.preStimAng = [Stats.cooled.FA.preStimAng;spikeField.precueMIFA.Angle(pos(n),n)];
        Stats.cooled.FA.postStimSPI = [Stats.cooled.FA.postStimSPI;max(spikeField.postcueMIFA.SPI(:,n))];
        Stats.cooled.FA.postStimAng = [Stats.cooled.FA.postStimAng;spikeField.postcueMIFA.Angle(pos(n),n)];
    end
end
end