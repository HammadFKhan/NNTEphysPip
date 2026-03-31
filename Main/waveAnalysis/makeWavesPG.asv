load("\\10.165.57.13\Sutter_backup\Hammad\Ephys\LeverTask\LeverTaskRebuttal\eOPN3\44263\Day7\44263_M_eOPN.mat");


% IntanBehaviour.MIHitTrace = [IntanBehaviourBaseline.MIHitTrace,IntanBehaviourOpto.MIHitTrace];
% IntanBehaviour.MIFATrace = [IntanBehaviourBaseline.MIFATrace,IntanBehaviourOpto.MIFATrace];


%% Making a new struct array to pass to the TCN 

% For CueHit conditions
WavesTCN.Hit = struct('levertrace', {IntanBehaviour.hitTrace.trace});
[WavesTCN.Hit.rawLFP] = deal(IntanBehaviour.hitTrace.rawLFP);

% For MIHit conditions
WavesTCN.MIHit = struct('levertrace', {IntanBehaviour.MIHitTrace.trace});
[WavesTCN.MIHit.rawLFP] = deal(IntanBehaviour.MIHitTrace.rawLFP);

% For MIFA conditions
WavesTCN.MIFA = struct('levertrace', {IntanBehaviour.MIFATrace.trace});
[WavesTCN.MIFA.rawLFP] = deal(IntanBehaviour.MIFATrace.rawLFP);


% Generating phase gradients 
filterOrder = 4;
filterLP1 = 5;
filterLP2 = 40;

% CueHit
for i=1:size(WavesTCN.MIHit,2)
    WavesTCN.Hit(i).xf = bandpass_filter(WavesTCN.Hit(i).rawLFP,filterLP1,filterLP2,filterOrder,parameters.Fs);
    [WavesTCN.Hit(i).xgp, WavesTCN.MIHit(i).wt] = generalized_phase(WavesTCN.Hit(i).xf,parameters.Fs,0);
end

% MIHit
for i=1:size(WavesTCN.MIHit,2)
    WavesTCN.MIHit(i).xf = bandpass_filter(WavesTCN.MIHit(i).rawLFP,filterLP1,filterLP2,filterOrder,parameters.Fs);
    [WavesTCN.MIHit(i).xgp, WavesTCN.MIHit(i).wt] = generalized_phase(WavesTCN.MIHit(i).xf,parameters.Fs,0);
end

% MIFA

for i=1:size(WavesTCN.MIFA,2)
    WavesTCN.MIFA(i).xf = bandpass_filter(WavesTCN.MIFA(i).rawLFP,filterLP1,filterLP2,filterOrder,parameters.Fs);
    [WavesTCN.MIFA(i).xgp, WavesTCN.MIFA(i).wt] = generalized_phase(WavesTCN.MIFA(i).xf,parameters.Fs,0);
end

% Getting PG and V
[WavesTCN.Hit] = getPG(WavesTCN.Hit,parameters);
[WavesTCN.MIHit] = getPG(WavesTCN.MIHit,parameters);
[WavesTCN.MIFA] = getPG(WavesTCN.MIFA,parameters);
WavesHit = WavesTCN.Hit;
%% Saving
savepath = uigetdir(path);
sessionName = [savepath,'/','NoTagSom_Day1GridsWaves.mat'];
save(sessionName,"WavesHit","IntanBehaviour","fpath","parameters","savepath","-v7.3");