function LFP = bestLFP(LFP)

pow = rms(LFP.LFP,2);
[~,loc] = max(pow);
disp(['Best LFP: ' num2str(loc)])
count = 1;
for chan = 1:size(LFP,1)
    rho = corr(LFP.LFP(1,:)',LFP.LFP(chan,:)');
    if rho>0.9
        medianLFP(count,:) = LFP.LFP(chan,:);
        count = count+1;
    end
end
disp('Resampling to original Frequency');
LFP.coherenceLFP = resample(median(medianLFP,1),8192,1024);


