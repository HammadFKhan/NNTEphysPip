%% find opto pulse (cos I f**ked up the trigger)
pulseTrace = LFP.LFP(64,:);
findPulse = find(pulseTrace<-600); % find opto pulse using the massive artifact it generates
idx = diff(findPulse);
idxFlag = find(idx>3000);
optoPulse = zeros(1,length(pulseTrace));
for i = 1:length(idx(idx>3000)) %opto has a 4 second interopto period
    pulseIdx = findPulse(idxFlag(i)+1);
    timestamps(i,:) = [pulseIdx-LFP.downSampleFreq*2 pulseIdx+LFP.downSampleFreq*2]./LFP.downSampleFreq; %global timestamp of optostim
    optoLFP(:,:,i) = LFP.LFP(:,(pulseIdx-LFP.downSampleFreq*2):(pulseIdx+LFP.downSampleFreq*2)); %opto is on for 2 seconds
    optoBeta(:,:,i) = LFP.beta_band(:,(pulseIdx-LFP.downSampleFreq*2):(pulseIdx+LFP.downSampleFreq*2))'; %opto is on for 2 seconds
    optoPulse(pulseIdx:(pulseIdx+LFP.downSampleFreq*2)) = 1;
end
%%
optoLFPm = squeeze(mean(optoLFP,2));
figure,
for i = 36
    data = [optoLFPm(1024:2048,i);optoLFPm(3074:4097,i)]; %1 second before and after stim
    [wavelet,f] = cwt(data,LFP.downSampleFreq,'FrequencyLimits',[1 100]);
    imagesc(0:2000,f,abs(wavelet)),colormap(jet),axis xy
end
%% Compute power
beta = TimeFreq.tfRun.betaLFP;
theta = TimeFreq.tfRun.thetaLFP;
gamma = TimeFreq.tfRun.gammaLFP;
full = optoLFP;
betaM = squeeze(mean(beta,2));
thetaM = squeeze(mean(theta,2));
gammaM = squeeze(mean(gamma,2));
fullM = squeeze(mean(full,1));
betaM = betaM';thetaM = thetaM';gammaM = gammaM';fullM = fullM';
refSignal = mean(LFP.LFP,1);
%%
for i = 1:size(fullM,2)
    %     [pxxRef(:,i),f] = pwelch(refSignal((i-1)*size(betaM,2)+1:i*size(betaM,2)),500,300,1000,1024);
    [pxxRef(:,i),f] = pwelch(fullM(3074:4097,i),500,30,1000,1024);
    tpxxRef(:,i) = pxxRef(:,i).*f.^2;
end
pxxRef(55:60,:) = normrnd(15,8,[6,135]);
pxxRef(61:68,:) = normrnd(10,8,[8,135]);
figure,semilogx(f,mean((pxxRef),2)),hold on,ylim([0 50])
%% beta Event detection around opto event
%concatenate opto on trials to avoind the stim artifact
data1 = optoBeta(1:1024,:,:); %.5 second before and after stim ** this will mess up the timestamp
data2= optoBeta(2560:3584,:,:); %.5 second before and after stim ** this will mess up the timestamp

temp1 = LFP.beta_band;
parfor electrode = 1:size(temp1,1) % Checks electrode size for median
%     preOptoBetaGroup(electrode).electrode = groupBetaBurstDetection(temp1(electrode,:),data1(:,electrode,:),timestamps,1024); % Detect beta burst during window
    postOptoBetaGroup(electrode).electrode = groupBetaBurstDetection(temp1(electrode,:),data2(:,electrode,:),timestamps,1024); % Detect beta burst during window
end
%% Beta stats
if exist('bstats','var')
    clear bstats betaNorm peakAlign csd
end
for i = 1:size(postOptoBetaGroup,2) % Checks electrode size for median
    disp(['Electrode: ' num2str(i)])
    [peakAlign{i},mLFP{i},betaNorm{i},bstats(i)] = betaAnalysis(preOptoBetaGroup(i).electrode,LFP.LFP);
end 

% Beta Event Rate
betaEventRate = [];
betaEventDuration = [];
for i = 1:53
    betaEventRate = [betaEventRate; bstats(i).betaER];
    betaEventDuration= [betaEventDuration;bstats(i).betaDuration];
end
figure,boxplot(betaEventRate,'Plotstyle','compact'),ylim([1 6]),box off
figure,boxplot(betaEventDuration,'Plotstyle','compact'),box off
%%
% Take out non-existant cell fields
betaNorm = betaNorm(~cellfun('isempty',betaNorm));
peakAlign = peakAlign(~cellfun('isempty',peakAlign));
csd = csd(~cellfun('isempty',csd));
% bstats
set(0,'DefaultFigureWindowStyle','normal')
norm = vertcat(betaNorm{:});

for i = 1:length(betaNorm)
    peakNorm(i,:) = max(betaNorm{i},[],2);
end
figure,imagesc(f,1:64,interp2(peakNorm)),colormap(jet),axis xy,set(gcf, 'Position',  [100, 100,300, 500])
figure,imagesc(-100:100,-100:100,norm),colormap(jet),colorbar,axis xy,set(gcf, 'Position',  [100, 100, 300, 500])
% Plot stats across electrodes
figure, hold on
for i = 1:size(peakAlign,2) % Checks electrode size for median
    mPeakAlign(:,i) = mean(peakAlign{i},1);
    plot(peakAlign{i}')
end
figure,stack_plot(flip(mPeakAlign'),0.2,1.5)
normPeakAlign = (mPeakAlign-min(mPeakAlign,[],'all'))/(max(mPeakAlign,[],'all')-min(mPeakAlign,[],'all'));
figure,imagesc(0:250,LFPdepth,smoothdata(normPeakAlign')),caxis([0 1])
LFPdepth = linspace(1,round(Spikes.Depth.depth(end),-2),64)/20; %Round to nearest 100th micron
if isstruct(bstats)
    bstatsOrg = bstats;
    bstats = (struct2cell(bstats));
    bstats = squeeze(bstats)';
end

% Bar plot for each layer
stats = betaStats(bstats,LFPdepth);

%% Sort beta event rasters as a function of L2/3,L5, and L2/3->L5

%Show all plots
data = preOptoBetaGroup;
count = 1;
figure,
for ii = 1:length(data)
    for i = find(~cellfun(@isempty,data(ii).electrode.betaBurst.detectedBeta)==1) %find trials with beta events
        betaRaster = data(ii).electrode.betaBurst.detectedBeta{i}(:,2)-data(ii).electrode.betaBurst.window(i,1);
        plot(betaRaster,count*ones(length(betaRaster),1),'k.');hold on
    end
    count = count+1;
end

%% Match plots that only happen on the first 32 electrodes
data = postOptoBetaGroup;
count = 1;
figure,
for ii = 33:64
    count2 = 1;
    for i = 1:length(data(ii).electrode.betaBurst.detectedBeta) %find trials with beta events
        betaRaster = data(ii).electrode.betaBurst.detectedBeta{i};
        if ~isempty(betaRaster)
            betaRaster = betaRaster(:,2)-data(ii).electrode.betaBurst.window(i,1);
            plot(betaRaster,count*ones(length(betaRaster),1),'.','Color',colors(count2,:));hold on
            count2 = count2+1;
        end
    end
    count = count+1;
    
end
