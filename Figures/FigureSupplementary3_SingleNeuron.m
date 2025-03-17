clear
%fpath = 'F:\LeverTask\Ephys\Analysis\M2Spikes';
%fpath = 'F:\LeverTask\Ephys\Analysis\spksPooledwFA'
%fpath = 'D:\M1_GSP';
fpath = 'D:\M2SpikeData';
file = dir(fullfile(fpath,'*.mat'));
for fileNum = 8:length(file)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(file(fileNum).folder,file(fileNum).name))
    Spikestotal(fileNum).fname = file(fileNum).name;
    Spikestotal(fileNum).Spikes = layerspikeAnalysis(Spikes,IntanBehaviour);
end
%%
hitFRStim = [];
hitFRnStim = [];
missFRStim = [];
missFRnStim = [];
hitProp = [];
missProp = [];
switchProp = [];
for n = 1:length(Spikestotal)
    dat = Spikestotal(n).Spikes.spikeProp;
    hitFRStim{n} = dat.hitFR(dat.hitspkTuning.stimResponsive,:);
    hitFRnStim{n} = dat.hitFR(~dat.hitspkTuning.stimResponsive,:);
    missFRStim{n} = dat.missFR(dat.missspkTuning.stimResponsive,:);
    missFRnStim{n} = dat.missFR(~dat.missspkTuning.stimResponsive,:);
    hitProp(n,1) = sum(dat.hitspkTuning.stimResponsive)/length(dat.hitspkTuning.stimResponsive);
    hitProp(n,2) = 1-hitProp(n,1);
    missProp(n,1) = sum(dat.missspkTuning.stimResponsive)/length(dat.missspkTuning.stimResponsive);
    missProp(n,2) = 1-missProp(n,1);
    switchProp(n,1) = sum(dat.hitspkTuning.stimResponsive==1 & dat.missspkTuning.stimResponsive==1)/length(dat.hitspkTuning.stimResponsive(dat.hitspkTuning.stimResponsive==1));
    switchProp(n,2) = 1-switchProp(n,1);
end

%%
cmap = cmocean('balance'); 
figure,subplot(221),imagesc(smoothdata(vertcat(hitFRStim{:}),2,'gaussian',20)),colormap(cmap),caxis([-3 3])
subplot(222),imagesc(smoothdata(vertcat(hitFRnStim{:}),2,'gaussian',20)),colormap(cmap),caxis([-3 3])
subplot(223),imagesc(smoothdata(vertcat(missFRStim{:}),2,'gaussian',20)),colormap(cmap),caxis([-3 3])
subplot(224),imagesc(smoothdata(vertcat(missFRnStim{:}),2,'gaussian',20)),colormap(cmap),caxis([-3 3])


figure,subplot(121),
pie(mean(hitProp,1))
title('Hit')
subplot(122),pie(mean(missProp,1))
title('Miss')
legend('Responsive','Unresponsive')
figure,
pie(nanmean(switchProp,1))
legend('Maintained Response','Became Unresponsive')