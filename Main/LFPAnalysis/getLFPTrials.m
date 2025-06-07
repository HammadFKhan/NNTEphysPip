function output = getLFPTrials(lfp,IntanBehaviour)
addpath(genpath('C:\Users\khan332\Documents\GitHub\generalized-phase'));
xo = lfp.LFP; %x,f1,f2,Fs
sz = size(xo);
xo = reshape(xo,sz(1),1,sz(2));

for i = 1:IntanBehaviour.nCueHit
    disp(['Analyzing trial: ' num2str(i)])
    hitWin = floor(IntanBehaviour.cueHitTrace(i).LFPtime*1000); %multiply by Fs;
    temp = xo(:,:,hitWin);
    output.hitLFP{i} = temp;
end

for i = 1:IntanBehaviour.nCueMiss
    disp(['Analyzing trial: ' num2str(i)])
    missWin = floor(IntanBehaviour.cueMissTrace(i).LFPtime*1000);
    temp = xo(:,:,missWin);
    output.missLFP{i} = temp;
end

%%% Motion Intiated
for i = 1:length(IntanBehaviour.MIHitTrace)
    disp(['Analyzing trial: ' num2str(i)])
    MIhitWin = floor(IntanBehaviour.MIHitTrace(i).LFPtime*1000); %multiply by Fs;
    temp = xo(:,:,MIhitWin);
    output.MIhitLFP{i} = temp;
end

for i = 1:length(IntanBehaviour.MIFATrace)
    disp(['Analyzing trial: ' num2str(i)])
    MIFAWin = floor(IntanBehaviour.MIFATrace(i).LFPtime*1000);
    temp = xo(:,:,MIFAWin);
    output.MIFALFP{i} = temp;
end
end