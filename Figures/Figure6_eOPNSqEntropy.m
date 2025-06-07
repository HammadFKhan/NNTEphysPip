ThalamuseOPN%% Sequentiality index for spiking activity during eopn perturbations
% Combining eOPN data together
M1eOPN = struct();
ThalamuseOPN = struct();

files = dir(fullfile('D:\eOPNData\M1Inactivation\','*.mat'));
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    M1eOPN(fileNum).IntanBehaviour = IntanBehaviour;
    M1eOPN(fileNum).filename = files(fileNum).name;
    [IntanBehaviourBaseline,IntanBehaviourOpto, Waves, WavesOpto] = separateOptoTrials(IntanBehaviour,IntanBehaviour.parameters);
    baselineId = 1:length(IntanBehaviourBaseline.cueHitTrace);
    eOPNId = length(IntanBehaviourBaseline.cueHitTrace)+1:length(IntanBehaviour.cueHitTrace);
    assert(length(eOPNId)==length(IntanBehaviourOpto.cueHitTrace))
    % Since the sqEntropy is accessing each spike structure we need to seperate the trials 
    baselineSpikes = Spikes; 
    eOPNSpikes = Spikes; 
    baselineSpikes.PSTH.hit.spks = cellfun(@(x) x(baselineId,:),baselineSpikes.PSTH.hit.spks,'UniformOutput',false);
    eOPNSpikes.PSTH.hit.spks = cellfun(@(x) x(eOPNId,:),eOPNSpikes.PSTH.hit.spks,'UniformOutput',false);

    M1eOPN(fileNum).baselineSqEntropy = getSqEntropy(baselineSpikes);
    M1eOPN(fileNum).eOPNSqEntropy = getSqEntropy(eOPNSpikes);
end

files = dir(fullfile('D:\eOPNData\ThalamusInactivation\','*.mat'));
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    ThalamuseOPN(fileNum).IntanBehaviour = IntanBehaviour;
    ThalamuseOPN(fileNum).filename = files(fileNum).name;

    [IntanBehaviourBaseline,IntanBehaviourOpto, Waves, WavesOpto] = separateOptoTrials(IntanBehaviour,IntanBehaviour.parameters);
    baselineId = 1:length(IntanBehaviourBaseline.cueHitTrace);
    eOPNId = length(IntanBehaviourBaseline.cueHitTrace)+1:length(IntanBehaviour.cueHitTrace);
    assert(length(eOPNId)==length(IntanBehaviourOpto.cueHitTrace))

    % Since the sqEntropy is accessing each spike structure we need to seperate the trials 
    baselineSpikes = Spikes; 
    eOPNSpikes = Spikes; 
    baselineSpikes.PSTH.hit.spks = cellfun(@(x) x(baselineId,:),baselineSpikes.PSTH.hit.spks,'UniformOutput',false);
    eOPNSpikes.PSTH.hit.spks = cellfun(@(x) x(eOPNId,:),eOPNSpikes.PSTH.hit.spks,'UniformOutput',false);

    ThalamuseOPN(fileNum).baselineSqEntropy = getSqEntropy(baselineSpikes);
    ThalamuseOPN(fileNum).eOPNSqEntropy = getSqEntropy(eOPNSpikes);
end

%%
dat = [];
for n = 1:length(M1eOPN)
    dat = vertcat(dat,[M1eOPN(n).baselineSqEntropy.CueHit.SqI(1)',M1eOPN(n).eOPNSqEntropy.CueHit.SqI(1)'])
end
figure,customBoxplot(dat)
[h,p] = ttest2(dat(:,1),dat(:,2))
ylabel('SI index')
title(['M1 pval = ',num2str(p)])

dat = [];
for n = 1:length(ThalamuseOPN)
    dat = vertcat(dat,[ThalamuseOPN(n).eOPNSqEntropy.CueHit.SqI(1)',ThalamuseOPN(n).baselineSqEntropy.CueHit.SqI(1)'])
end
figure,customBoxplot(dat)
[h,p] = ttest2(dat(:,1),dat(:,2))
ylabel('SI index')
title(['Thalamus pval = ',num2str(p)])