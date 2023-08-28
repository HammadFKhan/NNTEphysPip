%% Beta Event intracellular burst probability
isi = find(diff(trialFull(1,:))>4)/Fs;
isi = diff(isi);
isi = isi(isi>1/Fs);
%%
for nn = 1:size(trialFull,1)
    for n = 1:LFP(nn).trial.betaBurst.NumDetectedBeta
        isi = find(diff(trialFull(nn,:))>4)/Fs;
        temp = find(isi>=LFP(nn).trial.betaBurst.detectedBeta(n,1)  & isi<=LFP(nn).trial.betaBurst.detectedBeta(n,3));
        isi = diff(isi(temp));
        isi = isi(isi>(1/Fs));
        ISI(nn).burstProbability{n} = isi;
    end
    ISI(nn).betaAmplitude = mean(LFP(nn).trial.betaBurst.detectedBeta(:,4));
end
%%
burstProbability = [];
betaAmp = [];
for n = 1:size(trialFull,1)
    dat = ISI(n).burstProbability;
    len = length(dat);
    temp = cellfun(@(x) x<(1/50), dat, 'UniformOutput', false);
    temp = cellfun(@(x) sum(x),temp);
    temp = temp>0;
    temp = sum(temp)/len;
    burstProbability = vertcat(burstProbability,temp);
    betaAmp = vertcat(betaAmp,ISI(nn).betaAmplitude);
end
figure,histogram(burstProbability*100,0:5:100)