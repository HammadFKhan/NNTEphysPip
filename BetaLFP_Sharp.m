Fc = 3;
Wn = Fc./(Fs/2);
b1 = fir1(10000,Wn,'high');
disp('Filtering...')
tic;d = filtfilt(b1,1,lfp);toc;

Fc = [4 80];
Wn = Fc./(Fs/2);
b = fir1(10000,Wn,'bandpass');
tic;LFPfilt = filtfilt(b,1,d);toc;

filtered_data = customFilt(LFPfilt',Fs,[10 30]);
scaledLFP = filtered_data*1000;
LFP = IntrabetaBurstDetection(scaledLFP',Fs);
[peakAlignLFP,normLFP,fLFP,statsLFP] = IntrabetaAnalysis(LFP);

figure,plot(0:dt:(length(LFP.betaBurst.beta)-1)/Fs,LFP.betaBurst.beta);xlim([0 length(LFP.betaBurst.beta)/Fs]),box off
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

data = LFP2;
figure,
plot(0:dt:(length(data)-1)/Fs,data,'k');xlim([0 length(data)/Fs]),box off, hold on
for i = 1:size(LFP.betaBurst.detectedBeta,1)
xline(LFP.betaBurst.detectedBeta(i,2),'k');
end
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

