%% LFP-LFP coherance across trials between M2 and M1
% Here we only load in the data from simultenous M2 and M1 recordings which
% are pulled out from the CCA analysis
% Most of our data was only pertaining to spikes, we will instead make sure
% to build a pipeline grabbing LFP data only

rerun = 0; % We already have the data
if rerun
    % Load and parse matfile data
    % Combining eOPN data together
    lfpCoherence = struct();

    files = dir(fullfile('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\DualShank_LFP','*.mat'));
    for fileNum = 1:length(files)
        disp(['File number: ' num2str(fileNum)])
        load(fullfile(files(fileNum).folder,files(fileNum).name))
        lfpCoherence(fileNum).IntanBehaviour = IntanBehaviour;
        lfpCoherence(fileNum).filename = files(fileNum).name;

        xgp1 = LFP.probe1.genPhase.hitxgp; %% 64F P1
        xgp2 = LFP.probe2.genPhase.hitxgp; %% 64F P2
        xgp3 = LFP.probe3.genPhase.hitxgp; %% 64 Sharp P1
        [lfpCoherence(fileNum).phaseCor.preCuephaseCorrHit, lfpCoherence(fileNum).phaseCor.postCuephaseCorrHit] = getphaseLFP(xgp1,xgp3);

        xgp1 = LFP.probe1.genPhase.missxgp;
        xgp3 = LFP.probe3.genPhase.missxgp;
        [lfpCoherence(fileNum).phaseCor.preCuephaseCorrMiss, lfpCoherence(fileNum).phaseCor.postCuephaseCorrMiss] = getphaseLFP(xgp1,xgp3);

        xgp1 = LFP.probe1.genPhase.MIFAxgp;
        xgp3 = LFP.probe3.genPhase.MIFAxgp;
        [lfpCoherence(fileNum).phaseCor.preCuephaseCorrFA, lfpCoherence(fileNum).phaseCor.postCuephaseCorrFA] = getphaseLFP(xgp1,xgp3);
    end
    % Save data
    fpath = 'Y:\Hammad\Ephys\LeverTask\Data_for_Figures\DualShank_LFP\combined';
    sessionName = [fpath,'/','LFPCoherenceTotal.mat'];
    % save(sessionName,"IntanBehaviour","fpath","parameters","-v7.3");
    save(sessionName,"lfpCoherence","fpath","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
end
clear
load('Y:\Hammad\Ephys\LeverTask\Data_for_Figures\DualShank_LFP\combined\LFPCoherenceTotal.mat')
%% Plot out stats
% M2->M2: 1:32,1:32 electrodes
% M2->M1: 1:32,33:64 electrodes
% M1->M2: 33:64,1:32 electrodes
% M1->M1: 33:64,33:64 electrodes
M2M2 = [];
M2M1 = [];
M1M1 = [];
for n = 1:length(lfpCoherence)
    M2M2{n,1} = mean(lfpCoherence(n).phaseCor.preCuephaseCorrHit(1:32,1:32));
    M2M2{n,2} = mean(lfpCoherence(n).phaseCor.postCuephaseCorrHit(1:32,1:32));
    M2M1{n,1} = mean(lfpCoherence(n).phaseCor.preCuephaseCorrHit(1:32,33:64));
    M2M1{n,2} = mean(lfpCoherence(n).phaseCor.postCuephaseCorrHit(1:32,33:64));
    M1M1{n,1} = mean(lfpCoherence(n).phaseCor.preCuephaseCorrHit(33:64,33:64));
    M1M1{n,2} = mean(lfpCoherence(n).phaseCor.postCuephaseCorrHit(33:64,33:64));
end
%%
preCue = horzcat(M2M1{:,1});
postCue = horzcat(M2M1{:,2});
tot = [(horzcat(M1M1{:,2})-mean(horzcat(M1M1{:,1})))',(horzcat(M2M2{:,2})-mean(horzcat(M2M2{:,1})))',abs((postCue-mean(preCue)))'];
figure,
customBarplot(tot)
ylim([-1 1]),axis square
ylabel('Phase coherence')
set(gca,'tickdir','out','fontsize',12)
[p, tbl, stats] = anova1(tot, {'M1M1 Adjusted', 'M2M2 Adjusted', 'Abs PostCue Diff'}, 'off');
% Add ANOVA results to the plot
if p < 0.05
    % If ANOVA is significant, perform multiple comparison test
    [c,m,h,gnames] = multcompare(stats, 'Display','off');
    
    % Add ANOVA p-value to the top of the figure
    title_str = sprintf('Custom Bar Plot of Adjusted Data (ANOVA: p = %.4f)', p);
    title(title_str, 'FontSize', 8, 'FontWeight', 'bold');
end
%% Plot it out
preCuephaseCorr = lfpCoherence(1).phaseCor.preCuephaseCorrHit;
postCuephaseCorr = lfpCoherence(1).phaseCor.postCuephaseCorrHit;

load myMap
myMap = jet;
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

preCuephaseCorr = lfpCoherence(1).phaseCor.preCuephaseCorrMiss;
postCuephaseCorr = lfpCoherence(1).phaseCor.postCuephaseCorrMiss;

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

preCuephaseCorr = lfpCoherence(1).phaseCor.preCuephaseCorrFA;
postCuephaseCorr = lfpCoherence(1).phaseCor.postCuephaseCorrFA;

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
%%
lfpTrialsM2 = getLFPTrials(LFP.probe1,IntanBehaviour);
lfpTrialsM1 = getLFPTrials(LFP.probe3,IntanBehaviour);
%%
data1 = cellfun(@(x) squeeze(median(x,1)), lfpTrialsM1.missLFP, 'UniformOutput', false);
data1 = horzcat(data1{:});
data2 = cellfun(@(x) squeeze(median(x)), lfpTrialsM2.missLFP, 'UniformOutput', false);
data2 = horzcat(data2{:});

params.Fs = 1000;
params.tapers = [4 5];
params.fpass = [4 58];
params.trialave = 1;
[C,phi,S12,S1,S2,f]=coherencyc(data1,data2,params);
figure,plot(f,C)
%%
movingwin = [0.5 0.05];
[C,phi,S12,S1,S2,t,f]=cohgramc(data1,data2,movingwin,params);
figure,imagesc(mean(C,3)')
%%
figure,plotSpectrogram(C',t,f,'surf','Wavelet Based Spectrogram','Time (s)','Frequency (Hz)')
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

function plotSpectrogram(c,t,f,plottype,varargin)

params = parseinputs(varargin{:});

if strcmp(plottype,'contourf')
    contourf(t,f,c,40,'LineColor','none');
else
    surf(t,f,c,EdgeColor='none');
    view(0,90);
end
shading interp;colormap(jet);
axis tight;
h = colorbar;
h.Label.String = 'Power in dB';
if isempty(params.xlab) && isempty(params.ylab)
    ylabel('Frequency (Hz)');xlabel('Time (s)');
else
    xlabel(params.xlab); ylabel(params.ylab);
end
if ~isempty(params.PlotTitle)
    title(params.PlotTitle);
else
    title('Wavelet based Spectogram');
end

%----------------------------------------------------------------
function params = parseinputs(varargin)
    
    params.PlotTitle = [];
    params.xlab = [];
    params.ylab = [];

    if isempty(varargin)
        return;
    end
    Len = length(varargin);
    if (Len==1)
        params.PlotTitle = varargin{1};
    end
    if (Len == 3)
        params.PlotTitle = varargin{1};
        params.xlab = varargin{2};
        params.ylab = varargin{3};
    end
end
end
