% Code for cooling figure
%% For wave RT
load('Y:\Hammad\Ephys\LeverTask\Cooling\23960SomCooling\Day3\TWSpeed_RT.mat')
colors = [ [181 29 98]/255; [43 198 214]/255];
figure,violinplot(temp,[],'ShowData',true,'ShowWhiskers',false,'ShowBox',false,'MarkerSize',5,'ViolinColor' ,colors);
set(gca,'tickdir','out','fontsize',16),box off,ylabel('Wave Speed cm/s'),ylim([0 40])
%% Baseline data
PGD = arrayfun(@(x) x.PGD, Waves.wavesHit, 'UniformOutput', false);
PGD = vertcat(PGD{:});
waveTime = arrayfun(@(x) x.waveTime, Waves.wavesHit,'UniformOutput',false);
waveTime = vertcat(waveTime{:});
waveDir = arrayfun(@(x) x.waveDir, Waves.wavesHit,'UniformOutput',false);
waveDir = horzcat(waveDir{:})';
waveAmp = arrayfun(@(x) horzcat(x.waveAmp), Waves.wavesHit,'UniformOutput',false);
waveAmp = horzcat(waveAmp{:});
rhoWavePre = {};
rhoWavePost = {};
for n = 1:length(Waves.wavesHit)
    idx = Waves.wavesHit(n).evaluationPoints<1500;
    rhoWavePre{n} = cellfun(@mean,Waves.wavesHit(n).rho(idx));
    idx = Waves.wavesHit(n).evaluationPoints<2000 & Waves.wavesHit(n).evaluationPoints>1501;
    rhoWavePost{n} = cellfun(@mean,Waves.wavesHit(n).rho(idx));
end
rhoWave = [cellfun(@mean,rhoWavePre)' cellfun(@mean,rhoWavePost)'];
normrhoWave = [rhoWave(:,1)/mean(rhoWave(:,1))-1 rhoWave(:,2)/mean(rhoWave(:,1))-1]*100;
figure,
customBoxplot(normrhoWave),hold on
for n = 1:size(normrhoWave,1)
    line(1:2,[normrhoWave(n,1) normrhoWave(n,2)]);
end
box off,set(gca,'tickdir','out','fontsize',14)

PGDpreCue = mean(PGD(:,1:1500),2);
PGDpostCue = mean(PGD(:,1501:2000),2);
normPGDchange = abs([PGDpreCue/mean(PGDpreCue)-1 PGDpostCue/mean(PGDpreCue)-1])*100;
figure,
customBoxplot(normPGDchange),hold on
for n = 1:length(PGDpreCue)
    line(1:2,[normPGDchange(n,1) normPGDchange(n,2)]);
end
box off,set(gca,'tickdir','out','fontsize',14)
%%
[waveCool.waveDir,waveCool.waveRho,waveCool.rhoMax,waveCool.rhoDiff] = wavePresentRho(Waves);
%close all
%% Cooling data
PGDc = arrayfun(@(x) x.PGD, Wavesc.wavesHit, 'UniformOutput', false);
PGDc = vertcat(PGDc{:});
waveTimec = arrayfun(@(x) x.waveTime, Wavesc.wavesHit,'UniformOutput',false);
waveTimec = vertcat(waveTimec{:});
waveDirc = arrayfun(@(x) x.waveDir, Wavesc.wavesHit,'UniformOutput',false);
waveDirc = horzcat(waveDirc{:})';
waveAmpc = arrayfun(@(x) horzcat(x.waveAmp), Wavesc.wavesHit,'UniformOutput',false);
waveAmpc = horzcat(waveAmpc{:});
PGDcpreCue = mean(PGDc(:,1:1500),2);
PGDcpostCue = mean(PGDc(:,1501:2000),2);
normPGDcchange = [PGDcpreCue/mean(PGDcpreCue)-1 PGDcpostCue/mean(PGDcpreCue)-1]*100;
figure,
customBoxplot(normPGDcchange)
for n = 1:length(PGDcpreCue)
    line(1:2,[normPGDcchange(n,1) normPGDcchange(n,2)]);
end
box off,set(gca,'tickdir','out','fontsize',14)
%% calculate % change
t1 = normPGDchange(:,2)-normPGDchange(:,1);
t2 = normPGDcchange(:,2)-normPGDcchange(:,1);
%%
figure,histogram(baselineWaveTimeCue,0:2:80,'LineStyle','none','Normalization','probability'),hold on
histogram(cooledWaveTimeCue,0:2:80,'LineStyle','none','Normalization','probability')
box off,set(gca,'tickdir','out','fontsize',14),xlabel('Wave Duration (ms)'),ylabel('Frequency')
ranksum(baselineWaveTimeCue,cooledWaveTimeCue)