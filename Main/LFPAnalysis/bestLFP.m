function LFP = bestLFP(LFP)
tic
pow = rms(LFP.LFP,2);
[~,loc] = max(pow);
disp(['Best LFP: ' num2str(loc)])
count = 1;
for chan = 1:size(LFP.LFP,1)
    rho = corr(LFP.LFP(loc,:)',LFP.LFP(chan,:)');
    if rho>0.7
        medianLFP(count,:) = LFP.LFP(chan,:);
        count = count+1;
    end
end
disp('Resampling to original Frequency');
LFP.medianLFP = LFP.LFP;
LFP.bestLFP = median(medianLFP,1); %Best LFP for Phase Phase Analysis
commonModeAvg = medianLFP-mean(medianLFP);
LFP.commonModeAvg = commonModeAvg;
LFP.coherenceLFP = resample(LFP.bestLFP,8192,1024); %Spike coherence LFP
LFP.times = (1:size(LFP.bestLFP,2))/LFP.downSampleFreq;
toc



