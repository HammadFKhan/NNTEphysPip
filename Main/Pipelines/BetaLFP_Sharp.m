Fc = 3;
Wn = Fc./(Fs/2);
b1 = fir1(10000,Wn,'high');
disp('Filtering...')
tic;d = filtfilt(b1,1,lfp);toc;

Fc = [1 100];
Wn = Fc./(Fs/2);
b = fir1(10000,Wn,'bandpass');
tic;
LFPfilt = filtfilt(b,1,Intracellular);toc;

Fc = [3000];
Wn = Fc./(Fs/2);
b = fir1(20000,Wn,'low');
tic;
Intrafilt2 = filtfilt(b,1,Intracellular);toc;

LFPfilt = customFilt(lfp',Fs,[10 30]);
scaledLFP = LFPfilt*1000;
filtIntra = customFilt(Intracellular',Fs,[4 12]);
%% generate trial dependancy for thresholding sake
trialLFP = [];trialIntra = [];trialFull = [];
idx = 1:10:(length(scaledLFP))/Fs;
for i = 1:length(idx)-1
    trialLFP(i,:) = scaledLFP(idx(i)*Fs:idx(i+1)*Fs);
%     trialIntra(i,:) = Intrafilt(idx(i)*Fs:idx(i+1)*Fs);
    trialFull(i,:) = Intracellular(idx(i)*Fs:idx(i+1)*Fs);
end
%
%% LFP intra /extra coherence
params.Fs = Fs;
params.fpass = [1 100];
params.tapers = [5 9];
movingwin = [0.5 0.5];
params.pad = 1;
params.trialave = 1;
[C,phi,S12,S1,S2,f]=coherencyc(trialLFP',trialFull',params);
[S,f,Serr] = mtspectrumc(trialLFP',params);
figure,semilogy(f,S);hold on
semilogy(f,Serr)
figure,plot(f,smoothdata(mean(C,2)))
pxxRef = [];tpxxRef = [];
for i = 1:55
    %     [pxxRef(:,i),f] = pwelch(refSignal((i-1)*size(betaM,2)+1:i*size(betaM,2)),500,300,1000,1024);
    [pxxRef(:,i),f] = pwelch(trialIntra(i,:),500,100,10000,Fs);
    tpxxRef(:,i) = pxxRef(:,i).*f.^2;
end
% pxxRef(55:60,:) = normrnd(15,8,[6,135]);
% pxxRef(61:68,:) = normrnd(10,8,[8,135]);
figure,plot(f,mean(pxxRef,2)),hold on,xlim([0 50])

%%
for i = 1:size(trialLFP,1)
%     LFP(i).trial = IntrabetaBurstDetection(trialLFP(i,:)',Fs);
    LFP(i).trial = IntrabetaBurstDetection(trialLFP(i,:)',Fs);
end
%%
[peakAlignLFP,normLFP,fLFP,statsLFP] = IntrabetaAnalysis(LFP(1).trial);

figure,plot(0:dt:(length(LFP.betaBurst.beta)-1)/Fs,LFP.betaBurst.beta);xlim([0 length(LFP.betaBurst.beta)/Fs]),box off
figure,plot(0:dt:(length(Intracellular)-1)/Fs,Intracellular);xlim([0 length(LFP.betaBurst.beta)/Fs]),box off
figure,plot(0:dt:(length(Intra.betaBurst.beta)-1)/Fs,Intra.betaBurst.beta);xlim([0 length(Intra.betaBurst.beta)/Fs]),box off
%% Intracellular and LFP pipette recordings

% Extra/Intra beta event correlation

% Create a beta window
s1 = size(LFP.betaBurst.detectedBeta,1);
s2 = size(Intra.betaBurst.detectedBeta,1);

if s1>s2
    winTrials = s1;
    betaWin = LFP.betaBurst.detectedBeta;
else 
    winTrials = s2;
    betaWin = Intra.betaBurst.detectedBeta;
end
%%
count = 1;
timestamp = [];
lfpDurationStamp = [];
intraDurationStamp = [];
durationStamp = [];

for i = 1:winTrials
    win = [(betaWin(i,2)-.250) (betaWin(i,2)+.250)];
    idx = find(LFP.betaBurst.detectedBeta(:,2)>=win(1) & LFP.betaBurst.detectedBeta(:,2)<=win(2));
    if ~isempty(idx)
        timestamp{count} = LFP.betaBurst.detectedBeta(idx,2) - win(1);
        lfpDurationStamp{count} = LFP.betaBurst.detectedBeta(idx,3)-LFP.betaBurst.detectedBeta(idx,1);
        intraDurationStamp{count} = Intra.betaBurst.detectedBeta(i,3)-Intra.betaBurst.detectedBeta(i,1);
        durationStamp{count} = [mean(LFP.betaBurst.detectedBeta(idx,3)-LFP.betaBurst.detectedBeta(idx,1)) mean(Intra.betaBurst.detectedBeta(i,3)-Intra.betaBurst.detectedBeta(i,1))];
        count = count+1;
    end
end

%%
if iscell(timestamp) && iscell(intraDurationStamp) && iscell(lfpDurationStamp) && iscell(durationStamp)
    timestamp = vertcat(timestamp{:}); 
    intraDurationStamp = vertcat(intraDurationStamp{:});
    lfpDurationStamp = vertcat(lfpDurationStamp{:});
    durationStamp = vertcat(durationStamp{:});
end

figure,histogram(timestamp-.25,-.25:0.025:.25) % Plot correlation distribution to betaWin


figure,plot(Intra.betaBurst.detectedBeta(:,2),2*ones(261,1),'r.');s

%%


data = Intracellular2;
figure,
plot(0:dt:(length(data)-1)/Fs,data,'r');xlim([0 length(data)/Fs]),box off,hold on
for i = 1:size(Intra.betaBurst.detectedBeta,1)
xline(Intra.betaBurst.detectedBeta(i,2),'k');
end

data = Intra.betaBurst.beta.^2;
figure,
plot(0:dt:(length(data)-1)/Fs,data,'k');xlim([0 length(data)/Fs]),box off, hold on
for i = 1:size(LFP.betaBurst.detectedBeta,1)
xline(LFP.betaBurst.detectedBeta(i,2),'k');
end
%% Intra and extra cross-correlation
[c,lags] = xcorr(Intracellular2,Intra.betaBurst.beta,Fs/2);


% Normalize c from -1 to 1
cNorm = 2*((c-min(c))/(max(c)-min(c)))-1;
figure,plot(lags/Fs,cNorm),box off

[c1,lags1] = xcorr(Intracellular2,LFP.betaBurst.beta,Fs/2);
figure,plot(lags1/Fs,c1),box off

[c2,lags2] = xcorr(Intra.betaBurst.beta,LFP.betaBurst.beta,Fs/2);
c2Norm = 2*((c-min(c2))/(max(c2)-min(c2)))-1;
figure,plot(lags2/Fs,c2Norm),box off
%% Intracellular Beta Window
count = 1;
ISI = [];
for i = 1:size(Intra.betaBurst.detectedBeta,1)
    winAP = [floor((Intra.betaBurst.detectedBeta(i,2)-0.05)*Fs) floor((Intra.betaBurst.detectedBeta(i,2)+0.05)*Fs)];
    betaAP{i} = Intracellular2(winAP(1):winAP(2));
    % Calculate Intracellular ISI
    data = betaAP{i};
    data = data-min(data);
    thresh = find(data>2);
    tri = diff(thresh);
    trig = find(tri~=1);
    trig = [1;trig]; %pad for first index
    COM = thresh(trig+1,1);
    ISI = [ISI;diff(COM/Fs)*1000]; %ISI in milliseconds
    for ii = 1:length(trig)
        if (COM(ii)+(0.002*Fs))>length(data) || (COM(ii)-(0.002*Fs))<1
            count = count+1;
        else
            wint(:,count) = data(COM(ii)-(0.001*Fs):COM(ii)+(0.002*Fs),1);
            count = count+1;
        end
        
    end
    
end
%% Vm amplitude inside/outside beta event window
betaEventVmShuf = [];fullVm = {}; betaEventVm = {};
for i = 1:length(LFP)% trial
    for ii = 1:(LFP(i).trial.betaBurst.NumDetectedBeta-1)
        if ii<2
            continue;
        end
        win = [LFP(i).trial.betaBurst.detectedBeta(ii,2)-0.1, LFP(i).trial.betaBurst.detectedBeta(ii,2)+0.1];
        try
            betaEventVm{i}(ii,:) = trialIntra(i,(win(1))*Fs:((win(1))*Fs)+4000);
            betaEventWin{i}(ii,:) = trialLFP(i,(win(1))*Fs:((win(1))*Fs)+4000);
            fullVm{i}(ii,:) = trialFull(i,(win(1))*Fs:((win(1))*Fs)+4000);
        catch ME
            continue;
        end
    end
%     mBetaEventVm(i,:) = mean(betaEventVm,1);
%     figure,plot(mBetaEventVm(i,:)),hold on
end
%%
betaEventTotal = vertcat(betaEventWin{:});
figure,hold on
for i = 1:55
    subplot(8,7,i),plot(mean(betaEventVm{i},1)),title(num2str(i))
end

total = vertcat(betaEventVm{:});
total2 = (total-mean(total,2))*10;
figure,hold on
plot(mean(total,1)*10)

figure,lineError((0:4000)/Fs,total2)
%% under beta event window per trials (box plot)
%Intra and extra cross-correlation
mC = [];
for i = 1:size(betaEventTotal,1)
    [c,lags] = xcorr(total2(i,:),betaEventTotal(i,:),Fs/2,'normalize');
    mC(i,1) = mean(abs(c));
end

figure,boxplot(mC)
% Normalize c from -1 to 1
cNorm = 2*((c-min(c))/(max(c)-min(c)))-1;
figure,plot(lags/Fs,c),box off

[c1,lags1] = xcorr(VmIntraL23scaled,LFP(1).trial.betaBurst.beta,Fs/2);
figure,plot(lags1/Fs,c1),box off

[c2,lags2] = xcorr(Intra(1).betaBurst.beta,LFP(1).trial.betaBurst.beta,Fs/2);
figure,plot(lags2/Fs,c2),box off

%%
comVect = total(1:25,:);
[X1] = featureProject(comVect,5,0);
figure,
    scatter3(X1(:,1),X1(:,2),X1(:,3),10,[236 0 140]/255,'filled'); hold on; %[43 57 144]/255
%%
figure,histogram(l23_1,0:0.05:1,'normalization','probability')
hold on
histogram(l5_1,0:0.1:1,'normalization','probability')