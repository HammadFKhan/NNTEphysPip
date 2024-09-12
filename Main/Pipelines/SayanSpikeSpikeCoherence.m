% figure,hold on

data = [];
for i  = 1:length(spikeCoherence)
id = i;
% plotspikeCoherence(spikeCoherence(id).spikea,spikeCoherence(id).spikeb,spikeCoherence(id).spikecoherence)
data(i,:) = smoothdata(cell2mat(spikeCoherence(id).spikecoherence.coho),'gaussian',5);
depth(i,:) = [Spikes.Clusters(spikeCoherence(id).neuronA).spikeDepth,Spikes.Clusters(spikeCoherence(id).neuronB).spikeDepth];
end
f = spikeCoherence(1).spikecoherence.freq{1};
beta = data(:,f>=13 & f<30);
gamma = data(:,f>=30 & f<90);
theta = data(:,f>=5 & f<12);
% %%
% figure,subplot(131),plot(theta','k'), hold on
% subplot(132),plot(beta','r'),hold on
% subplot(133),plot(gamma','b'),hold on
% %%
% figure,
% customBoxplot([mean(theta,2) mean(beta,2) mean(gamma,2)])
%% curate data
if ~exist('totCoherence','var'), totCoherence = [];end
if ~exist('totDepth','var'), totDepth = [];end

totCoherence = vertcat(totCoherence,data);
totDepth = vertcat(totDepth,depth);
% figure,plot(mean(data(id,:),1))
% %% 
temp = mean(SpikeCoherenceStats.Mon.coherence.L23L5beta,2);
id = find(temp>=0.012);
tot = vertcat(tot,data(id,:));
%% do stats
Monbeta = Montot(:,f>=13 & f<30);
Mongamma = Montot(:,f>=30 & f<90);
Montheta = Montot(:,f>=5 & f<12);
CTXbeta = CTXtot(:,f>=13 & f<30);
CTXgamma = CTXtot(:,f>=30 & f<90);
CTXtheta = CTXtot(:,f>=5 & f<12);
STRbeta = STRtot(:,f>=13 & f<30);
STRgamma = STRtot(:,f>=30 & f<90);
STRtheta = STRtot(:,f>=5 & f<12);
%% do stats with layers
Monbeta = MontotCoherence(:,f>=13 & f<30);
Mongamma = MontotCoherence(:,f>=30 & f<90);
Montheta = MontotCoherence(:,f>=5 & f<12);
CTXbeta = CTXtotCoherence(:,f>=13 & f<30);
CTXgamma = CTXtotCoherence(:,f>=30 & f<90);
CTXtheta = CTXtotCoherence(:,f>=5 & f<12);
STRbeta = STRtotCoherence(:,f>=13 & f<30);
STRgamma = STRtotCoherence(:,f>=30 & f<90);
STRtheta = STRtotCoherence(:,f>=5 & f<12);
%% now do by depth (L23-L23 L23-L5 L5-L5)
SpikeCoherenceStats.Mon.coherence = freqlayerCoherence(Montheta,Monbeta,Mongamma,MontotDepth);
SpikeCoherenceStats.CTX.coherence = freqlayerCoherence(CTXtheta,CTXbeta,CTXgamma,CTXtotDepth*25);
SpikeCoherenceStats.STR.coherence = freqlayerCoherence(STRtheta,STRbeta,STRgamma,STRtotDepth*25);

%%
figure,customBoxplot([mean(Montheta,2) mean(CTXtheta,2)]),title('CTX Theta')
figure,customBoxplot([mean(Monbeta,2) mean(CTXbeta,2)]),title('CTX Beta')
figure,customBoxplot([mean(Mongamma,2) mean(CTXgamma,2)]),title('CTX Gamma')
figure,customBoxplot([mean(Montheta,2) mean(STRtheta,2)]),title('STR Theta')
figure,customBoxplot([mean(Monbeta,2) mean(STRbeta,2)]),title('STR Beta')
figure,customBoxplot([mean(Mongamma,2) mean(STRgamma,2)]),title('STR Gamma')
%%
figure,customBoxplot([mean(Montheta,2) mean(Monbeta,2) mean(Mongamma,2)]),title('Mon'),ylim([0 0.2]),box off,set(gca,'fontsize',16),set(gca,'tickdir','out')
figure,customBoxplot([mean(CTXtheta,2) mean(CTXbeta,2) mean(CTXgamma,2)]),title('CTX'),ylim([0 0.2]),box off,set(gca,'fontsize',16),set(gca,'tickdir','out')
figure,customBoxplot([mean(STRtheta,2) mean(STRbeta,2) mean(STRgamma,2)]),title('STR'),ylim([0 0.2]),box off,set(gca,'fontsize',16),set(gca,'tickdir','out')
%% Layer specific
plotfreqlayerCoherence(SpikeCoherenceStats.Mon.coherence)
%%
plotfreqlayerCoherence(SpikeCoherenceStats.STR.coherence)
%%
plotfreqlayerCoherence(SpikeCoherenceStats.CTX.coherence)
%% L23
figure,customBoxplot([mean(MonL23Theta,2) mean(MonL23Beta,2) mean(MonL23Gamma,2)]),title('Mon'),ylim([0 0.2]),box off,set(gca,'fontsize',16),set(gca,'tickdir','out')
figure,customBoxplot([mean(CTXL23Theta,2) mean(CTXL23Beta,2) mean(CTXL23Gamma,2)]),title('CTX'),ylim([0 0.2]),box off,set(gca,'fontsize',16),set(gca,'tickdir','out')
figure,customBoxplot([mean(STRL23Theta,2) mean(STRL23Beta,2) mean(STRL23Gamma,2)]),title('STR'),ylim([0 0.2]),box off,set(gca,'fontsize',16),set(gca,'tickdir','out')
figure,customBoxplot([mean(MonL23ThetaCTX,2) mean(MonL23BetaCTX,2) mean(MonL23GammaCTX,2)]),title('Mon'),ylim([0 0.2]),box off,set(gca,'fontsize',16),set(gca,'tickdir','out')

%% local functions
function plotfreqlayerCoherence(stats)
figure,subplot(3,1,[1 3]),customBoxplot([max(stats.L23L23theta,[],2)...
    max(stats.L23L23beta,[],2) max(stats.L23L23gamma,[],2)]),box off,set(gca,'Tickdir','out'),title('L23<->L23')
% subplot(312),customBoxplot([max(stats.L23L5theta,[],2)...
%     max(stats.L23L5beta,[],2) max(stats.L23L5gamma,[],2)]),box off,set(gca,'Tickdir','out'),title('L23<->L5')
% subplot(313),customBoxplot([max(stats.L5L5theta,[],2)...
%     max(stats.L5L5beta,[],2) max(stats.L5L5gamma,[],2)]),box off,set(gca,'Tickdir','out'),title('L5<->L5')
end
function output = freqlayerCoherence(thetacoherence,betacoherence,gammacoherence,depth)
idx = find(depth(:,1)<400 & depth(:,2)<400);
L23L23theta = thetacoherence(idx,:);
L23L23beta  = betacoherence(idx,:);
L23L23gamma = gammacoherence(idx,:);
idx = find((depth(:,1)<400 & depth(:,2)>400) | (depth(:,1)>400 & depth(:,2)<400));
L23L5theta = thetacoherence(idx,:);
L23L5beta = betacoherence(idx,:);
L23L5gamma = gammacoherence(idx,:);
idx = find(depth(:,1)>=400 & depth(:,2));
L5L5theta = thetacoherence(idx,:);
L5L5beta = betacoherence(idx,:);
L5L5gamma = gammacoherence(idx,:);

output.L23L23theta = L23L23theta;
output.L23L23beta = L23L23beta;
output.L23L23gamma = L23L23gamma;
output.L23L5theta = L23L5theta;
output.L23L5beta = L23L5beta;
output.L23L5gamma = L23L5gamma;
output.L5L5theta = L5L5theta;
output.L5L5beta = L5L5beta;
output.L5L5gamma = L5L5gamma;
end

function plotspikeCoherence(spikea,spikeb,spikeCoherence)

interval = 0:2000;
spiker = spikea;

figure(1);hold on
%****************
disp('Plotting rastergrams (slow) ...');
subplot('position',[0.1 0.4 0.4 0.55]);
make_nice_spike_raster(spiker);
V = axis;
axis([0 size(spikea{1},2) V(3) V(4)]);
grid on;
ylabel('Trial Number');
title(fprintf('Unit rasters'));
%*************
subplot('position',[0.1 0.1 0.4 0.3]);
smooth_window = 25;  % give sigma of 12.5ms
make_nice_mean_raster(spiker,smooth_window);
V = axis;
axis([0 size(spikea{1},2) V(3) V(4)]);
plot([interval(1),interval(1)],[V(3),V(4)],'k-'); hold on;
plot([interval(end),interval(end)],[V(3),V(4)],'k-'); hold on;
ylabel('Firing Rate');
xlabel('Time (ms)');

spiker = spikeb;
%****************
disp('Plotting rasters2 per trial (slow) ...');
subplot('position',[0.57 0.4 0.4 0.55]);
make_nice_spike_raster(spiker);
V = axis;
axis([0 size(spikeb{1},2) V(3) V(4)]);
grid on;
ylabel('Trial Number');
title(fprintf('Unit rasters'));
%*************
subplot('position',[0.57 0.1 0.4 0.3]);
smooth_window = 25;  % give sigma of 12.5ms
make_nice_mean_raster(spiker,smooth_window);
V = axis;
axis([0 size(spikeb{1},2) V(3) V(4)]);
plot([interval(1),interval(1)],[V(3),V(4)],'k-'); hold on;
plot([interval(end),interval(end)],[V(3),V(4)],'k-'); hold on;
ylabel('Firing Rate');
xlabel('Time (ms)');

 figure(2);hold on
colo = 'rbgy';
CNUM=1;
subplot(2,2,1);
for ii = 1:CNUM
    H = semilogx(spikeCoherence.freq{ii},spikeCoherence.coho{ii},[colo(ii),'-']); hold on;
    set(H,'Linewidth',2);
    H = semilogx(spikeCoherence.freq{ii},(spikeCoherence.coho{ii}+spikeCoherence.scoho{ii}),[colo(ii),':']); hold on;
    H = semilogx(spikeCoherence.freq{ii},(spikeCoherence.coho{ii}-spikeCoherence.scoho{ii}),[colo(ii),':']);
    H = semilogx(spikeCoherence.freq{ii},spikeCoherence.rcoho{ii},[colo(ii),'--']); hold on;
    set(H,'Linewidth',1);
end
ylabel('Coherence Magnitude');
xlabel('Frequency (Hz)');
title(sprintf('Spike A with Spike B'));

subplot(2,2,2);
for ii = 1:CNUM
    H = plot(spikeCoherence.freq{ii},spikeCoherence.phaso{ii},[colo(ii),'-']); hold on;
    set(H,'Linewidth',2);
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.phaso{ii}+spikeCoherence.sphaso{ii}),[colo(ii),':']); hold on;
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.phaso{ii}-spikeCoherence.sphaso{ii}),[colo(ii),':']);
    %H = plot(freq{ii},rphaso{ii},[colo(ii),'--']); hold on;
    set(H,'Linewidth',1);
end
ylabel('Coherence Angle');
xlabel('Freq');
title(sprintf('Spike A with Spike B'));

subplot(2,2,3);
for ii = 1:CNUM
    H = plot(spikeCoherence.freq{ii},spikeCoherence.spika_pow{ii},[colo(ii),'-']); hold on;
    set(H,'Linewidth',2);
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.spika_pow{ii}+spikeCoherence.sspika_pow{ii}),[colo(ii),':']); hold on;
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.spika_pow{ii}-spikeCoherence.sspika_pow{ii}),[colo(ii),':']);
end
ylabel('Spike Unit A Pow');
xlabel('Frequency (Hz)');
title(sprintf('Spike A with Spike B'));

subplot(2,2,4);
for ii = 1:CNUM
    H = plot(spikeCoherence.freq{ii},spikeCoherence.spike_pow{ii},[colo(ii),'-']); hold on;
    set(H,'Linewidth',2);
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.spike_pow{ii}+spikeCoherence.sspike_pow{ii}),[colo(ii),':']); hold on;
    H = plot(spikeCoherence.freq{ii},(spikeCoherence.spike_pow{ii}-spikeCoherence.sspike_pow{ii}),[colo(ii),':']);
end
ylabel('Spike Unit B Pow');
xlabel('Frequency (Hz)');
title(sprintf('Spike B with Spike A'));

end
%****************** function to make a nice raster plot ****************
function make_nice_mean_raster(spmat,smooth_window)
%*********** spmat1 and spmat2 are spike matrices of two conditions you wish to compare
%*********** smooth_window ... gaussian smoothing in millisecs
numconds = size(spmat,2);
if (numconds==2)
    colo = [[1,0,0];[0,0,1]];
else
    colo = [[1,0,0];[0,0,1];[0,1,0];[1,1,0]];
end

for k = 1:numconds
    spud = spmat{1,k};
    numtrials(k) = size(spud,1);
    smorate = gauss_smooth(sum( spud(1:numtrials(k),:))/....
        numtrials(k),smooth_window)*1000;
    H = plot(smorate,'k-'); hold on;
    set(H,'Color',colo(k,:));
end

end


%******************* make a rastergram of the actual spikes on each trial
function make_nice_spike_raster(spmat)

colo = [[1,1,1];[0,0,0]];
colormap(colo);

totspike = [];
for cubo = 1:size(spmat,2)
    totspike = [totspike; (cubo*spmat{1,cubo}) ];
    totspike = [totspike; (cubo*spmat{1,cubo}) ];
    
end
totspike = totspike + 1;
imagesc(totspike);
V = axis;
axis([V(1) V(2) 0 size(totspike,1)]);

end

%**************************************************************
function output = gauss_smooth(input, window)
% Smoothing function:
% output = smooth(input, window)
% "Window" is the total kernel width.
% Input array must be one-dimensional.

input_dims = ndims(input);
input_size = size(input);
if input_dims > 2 | min(input_size) > 1,
    disp('Input array is too large.');
    return
end

if input_size(2) > input_size(1),
    input = input';
    toggle_dims = 1;
else
    toggle_dims = 0;
end

if window/2 ~= round(window/2),
    window = window + 1;
end
halfwin = window/2;

input_length = length(input);
%********* gauss window +/- 1 sigma
x = -halfwin:1:halfwin;
kernel = exp(-x.^2/(window/2)^2);
kernel = kernel/sum(kernel);

padded(halfwin+1:input_length+halfwin) = input;
padded(1:halfwin) = ones(halfwin, 1)*input(1);
padded(length(padded)+1:length(padded)+halfwin) = ones(halfwin, 1)*input(input_length);

output = conv(padded, kernel);
output = output(window:input_length+window-1);

if toggle_dims == 1,
    output = output';
end
end