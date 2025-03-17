%% Combined wave dynamics with stats for figure
% Load data
if ~exist('wavesHit','var')
    load('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\M1GridBaseline\WavePropCombined.mat')
end

% Plot wave present scatter plot and histogram for trial types
datHit = []; datMiss = []; datMIHit = [];datFA = [];
dat = arrayfun(@(x) vertcat(x.wavePresent), wavesHit, 'UniformOutput', false);
datHit = vertcat(dat{:});
dat = arrayfun(@(x) vertcat(x.wavePresent), wavesMiss, 'UniformOutput', false);
datMiss = vertcat(dat{:});
dat = arrayfun(@(x) vertcat(x.wavePresent), wavesMIHit, 'UniformOutput', false);
datMIHit = vertcat(dat{:});
dat = arrayfun(@(x) vertcat(x.wavePresent), wavesMIFA, 'UniformOutput', false);
datFA = vertcat(dat{:});
%% plot rastergram(the data is too large here tho)
colors = [0 0.4470 0.7410;0.75 0.75 0.75;0 0.4470 0.7410;190/255 30/255 45/255]; % data is indexed as Hit, Miss, MIHIT, MIFA
count = 1;
time = (-1500:1500)/1000;
figure(1),clf
bin = 10;
for n = 1:bin:size(datHit,1)
    plot(time,count*datHit(n,:),'.','Color',colors(1,:)); hold on
    count = count+1;
end
for n = 1:bin:size(datMiss,1)
    plot(time,count*datMiss(n,:),'.','Color',colors(2,:)); hold on
    count = count+1;
end
for n = 1:bin:size(datMIHit,1)
    plot(time,count*datMIHit(n,:),'.','Color',colors(3,:)); hold on
    count = count+1;
end
for n = 1:bin:size(datFA,1)
    plot(time,count*datFA(n,:),'.','Color',colors(4,:)); hold on
    count = count+1;
end
%%
figure(2),clf
dat = smoothdata(datHit,2,'gaussian',5);
subplot(311),
bar(time,sum(dat,1)/size(dat,1),'EdgeColor','none','FaceColor',[0.8 0.8 0.8]), hold on
plot(time,smoothdata(sum(dat,1)/size(dat,1),'gaussian',100),'Color',colors(1,:),'LineWidth',2)
set(gca,'tickdir','out','fontsize',14),box off
xlim([-0.5 1.5])
dat = smoothdata(datMiss,2,'gaussian',5);
subplot(312),
bar(time,sum(dat,1)/size(dat,1),'EdgeColor','none','FaceColor',[0.8 0.8 0.8]), hold on
plot(time,smoothdata(sum(dat,1)/size(dat,1),'gaussian',50),'Color',colors(2,:),'LineWidth',2)
set(gca,'tickdir','out','fontsize',14),box off
xlim([-0.5 1.5])

dat = smoothdata(datFA,2,'gaussian',5);
subplot(313),
bar(time,sum(dat,1)/size(dat,1),'EdgeColor','none','FaceColor',[0.8 0.8 0.8]), hold on
plot(time,smoothdata(sum(dat,1)/size(dat,1),'gaussian',50),'Color',colors(4,:),'LineWidth',2)
set(gca,'tickdir','out','fontsize',14),box off
xlim([-0.5 1.5])
%%
dat = arrayfun(@(x) vertcat(x.PGD), wavesHit, 'UniformOutput', false);
dat = vertcat(dat{:});
% Stim normalize
dat  = mean(dat);
[val,id] = min(dat(1:1500));
dat(1:id-1) = dat(1:id-1)-((max(dat(:,1:1500))-val))/1.5;
figure(3),clf
plot(time,smoothdata(mean(dat,1),'gaussian',25))
set(gca,'tickdir','out','fontsize',14),box off
xlim([-0.5 1.5]),hold on

dat = arrayfun(@(x) vertcat(x.PGD), wavesMiss, 'UniformOutput', false);
dat = vertcat(dat{:});
% Stim normalize
dat  = mean(dat);
[val,id] = min(dat(1:1500));
dat(1:id-1) = dat(1:id-1)-((max(dat(:,1:1500))-val))/1.5;
plot(time,smoothdata(mean(dat,1),'gaussian',25))
set(gca,'tickdir','out','fontsize',14),box off
xlim([-0.5 1.5])

dat = arrayfun(@(x) vertcat(x.PGD), wavesMIFA, 'UniformOutput', false);
dat = vertcat(dat{:});
% Stim normalize
dat  = mean(dat);
[val,id] = min(dat(1:1500));
%dat(1:id-1) = dat(1:id-1)-((max(dat(:,1:1500))-val))/1.5;
figure(5),clf
plot(time,smoothdata(mean(dat,1),'gaussian',25))
set(gca,'tickdir','out','fontsize',14),box off
xlim([-0.5 1.5])

dat = arrayfun(@(x) vertcat(x.PGD), wavesFA, 'UniformOutput', false);
dat = vertcat(dat{:});
% Stim normalize
dat  = mean(dat);
[val,id] = min(dat(1:1500));
%dat(1:id-1) = dat(1:id-1)-((max(dat(:,1:1500))-val))/1.5;
figure(6),clf
plot(time,smoothdata(mean(dat,1),'gaussian',25))
set(gca,'tickdir','out','fontsize',14),box off
xlim([-0.5 1.5])
