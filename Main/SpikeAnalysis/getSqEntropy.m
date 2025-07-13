function SqEntropy = getSqEntropy(Spikes)
%% % Loop across bins to test temporal sparsity response (ie. we care about
% peak entropy which is congruent to precision based on binning)
% Calculate sequencial entropy as this is useful for check MIHIT to CueHit
% data which is not captured in the initial figure 1. 
CueHit = struct(); CueMiss = struct(); MIFA = struct(); MIHit = struct();
bin = 50:200:1500;
for n = 1:length(bin)
    disp(['Calculating with bin length: ' num2str(bin(n))])
    NumEntropyBins = bin(n);
    dat = Spikes.PSTH.hit.spks;
    Data = cell2mat(permute(dat,[1,3,2]));
    Data = permute(Data,[1,3,2]);
    [CueHit.SqI(n), CueHit.PE(n), CueHit.TS(n)] = SeqIndexDB(Data,NumEntropyBins);
    CueHit.bin = bin;
    
    dat = Spikes.PSTH.miss.spks;
    Data = cell2mat(permute(dat,[1,3,2]));
    Data = permute(Data,[1,3,2]);
    [CueMiss.SqI(n), CueMiss.PE(n), CueMiss.TS(n)] = SeqIndexDB(Data,NumEntropyBins);
    CueMiss.bin = bin;
    
    dat = Spikes.PSTH.MIFA.spks;
    Data = cell2mat(permute(dat,[1,3,2]));
    Data = permute(Data,[1,3,2]);
    [MIFA.SqI(n),MIFA.PE(n), MIFA.TS(n)] = SeqIndexDB(Data,NumEntropyBins);
    MIFA.bin = bin;
    
    dat = Spikes.PSTH.MIHit.spks;
    Data = cell2mat(permute(dat,[1,3,2]));
    Data = permute(Data,[1,3,2]);
    [MIHit.SqI(n), MIHit.PE(n), MIHit.TS(n)] = SeqIndexDB(Data,NumEntropyBins);
    MIHit.bin = bin;
end

SqEntropy.CueHit = CueHit;
SqEntropy.CueMiss = CueMiss;
SqEntropy.MIFA = MIFA;
SqEntropy.MIHit = MIHit;
%% Plot it out

% figure;
% subplot(131)
% plot(CueHit.bin,CueHit.SqI),hold on
% plot(CueHit.bin,CueMiss.SqI)
% plot(CueHit.bin,MIFA.SqI)
% plot(CueHit.bin,MIHit.SqI)
% axis square
% 
% subplot(132)
% plot(CueHit.bin,CueHit.PE),hold on
% plot(CueHit.bin,CueMiss.PE)
% plot(CueHit.bin,MIFA.PE)
% plot(CueHit.bin,MIHit.PE)
% axis square
% 
% subplot(133)
% plot(CueHit.bin,CueHit.TS),hold on
% plot(CueHit.bin,CueMiss.TS)
% plot(CueHit.bin,MIFA.TS)
% plot(CueHit.bin,MIHit.TS)
% axis square
