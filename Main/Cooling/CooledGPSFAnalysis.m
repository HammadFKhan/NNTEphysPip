function Spikes = CooledGPSFAnalysis(Spikes,lfp,IntanBehaviour)

% Assertion check
assert(size(Spikes.PSTH.hit.spks{1},2)==size(squeeze(lfp.hitxgp{1}),2)); %checks length
assert(size(Spikes.PSTH.hit.spks{1},1)==length(lfp.hitxgp));
assert(size(Spikes.PSTH.miss.spks{1},1)==length(lfp.missxgp));
assert(size(Spikes.PSTH.MIHit.spks{1},1)==length(lfp.MIhitxgp));
assert(size(Spikes.PSTH.MIFA.spks{1},1)==length(lfp.MIFAxgp));

% Define cooling trials 
lfp.Cooled.hitxgp = lfp.hitxgp(IntanBehaviour.hitTemp<-10);
lfp.Cooled.missxgp = lfp.missxgp(IntanBehaviour.missTemp<-10);
lfp.Cooled.MIhitxgp = lfp.MIhitxgp(IntanBehaviour.hitTemp<-10);
lfp.Cooled.MIFAxgp = lfp.MIFAxgp(IntanBehaviour.FATemp<-10);

lfp.Baseline.hitxgp = lfp.hitxgp(IntanBehaviour.hitTemp>-10);
lfp.Baseline.missxgp = lfp.missxgp(IntanBehaviour.missTemp>-10);
lfp.Baseline.MIhitxgp = lfp.MIhitxgp(IntanBehaviour.hitTemp>-10);
lfp.Baseline.MIFAxgp = lfp.MIFAxgp(IntanBehaviour.FATemp>-10);

Spikes.Cooled.PSTH.hit.spks = cellfun(@(x) x(IntanBehaviour.hitTemp<-10,:),Spikes.PSTH.hit.spks, 'UniformOutput', false);
Spikes.Cooled.PSTH.miss.spks = cellfun(@(x) x(IntanBehaviour.missTemp<-10,:),Spikes.PSTH.miss.spks, 'UniformOutput', false);
Spikes.Cooled.PSTH.MIHit.spks = cellfun(@(x) x(IntanBehaviour.hitTemp<-10,:),Spikes.PSTH.MIHit.spks, 'UniformOutput', false);
Spikes.Cooled.PSTH.MIFA.spks = cellfun(@(x) x(IntanBehaviour.FATemp<-10,:),Spikes.PSTH.MIFA.spks, 'UniformOutput', false);

Spikes.Baseline.PSTH.hit.spks = cellfun(@(x) x(IntanBehaviour.hitTemp>-10,:),Spikes.PSTH.hit.spks, 'UniformOutput', false);
Spikes.Baseline.PSTH.miss.spks = cellfun(@(x) x(IntanBehaviour.missTemp>-10,:),Spikes.PSTH.miss.spks, 'UniformOutput', false);
Spikes.Baseline.PSTH.MIHit.spks = cellfun(@(x) x(IntanBehaviour.hitTemp>-10,:),Spikes.PSTH.MIHit.spks, 'UniformOutput', false);
Spikes.Baseline.PSTH.MIFA.spks = cellfun(@(x) x(IntanBehaviour.FATemp>-10,:),Spikes.PSTH.MIFA.spks, 'UniformOutput', false);

Spikes.Baseline.Clusters = Spikes.Clusters;
Spikes.Cooled.Clusters = Spikes.Clusters;

% Then we just call the function twice here, jjeje
Spikes.Baseline = GPSFAnalysis(Spikes.Baseline,lfp.Baseline);
Spikes.Cooled = GPSFAnalysis(Spikes.Cooled,lfp.Cooled);
end
