function [IntanBehaviour,IntanBehaviourOpto, Waves, WavesOpto] = separateOptoTrials(IntanBehaviour,parameters,Waves)

if ~exist('Waves','var')
    wavesPassed = 0;
    disp('No waves passed.')
    Waves = [];
    WavesOpto = [];
else
    wavesPassed = 1;
end

%% Separate trials based on the opto 

firstOptoIndex = find(IntanBehaviour.optoTrace == 1,1);

if parameters.opto == 1
    % Cue Hit and MIHit
    firstOptoTrial = find(IntanBehaviour.cueHit(:,1)>firstOptoIndex,1);
    a = num2cell(zeros(1,firstOptoTrial));
    [IntanBehaviour.cueHitTrace(1:firstOptoTrial).opto] = a{:};
    [IntanBehaviour.MIHitTrace(1:firstOptoTrial:end).opto] = a{:};
    [IntanBehaviour.hitTrace(1:firstOptoTrial:end).opto] = a{:};

    a = num2cell(ones(1,IntanBehaviour.nCueHit-firstOptoTrial+1));
    [IntanBehaviour.cueHitTrace(firstOptoTrial:end).opto] = a{:};
    [IntanBehaviour.MIHitTrace(firstOptoTrial:end).opto] = a{:};
    [IntanBehaviour.hitTrace(firstOptoTrial:end).opto] = a{:};
    
    % Cue Miss
    firstOptoTrial = find(IntanBehaviour.cueMiss(:,1)>firstOptoIndex,1);
    a = num2cell(ones(1,IntanBehaviour.nCueMiss-firstOptoTrial+1));
    [IntanBehaviour.cueMissTrace(firstOptoTrial:end).opto] = a{:};
    a = num2cell(zeros(1,firstOptoTrial));
    [IntanBehaviour.cueMissTrace(1:firstOptoTrial).opto] = a{:};

    % FA and MIFA
    MIFAIndex = cell2mat(arrayfun(@(s) s.LFPIndex(parameters.windowBeforeMI*parameters.Fs+1),IntanBehaviour.MIFATrace,'UniformOutput',false))';
    firstOptoTrial = find(MIFAIndex>firstOptoIndex,1);
    a = num2cell(ones(1,IntanBehaviour.nMiss-firstOptoTrial+1));
    [IntanBehaviour.missTrace(firstOptoTrial:end).opto] = a{:};
    [IntanBehaviour.MIFATrace(firstOptoTrial:end).opto] = a{:};
    a = num2cell(zeros(1,firstOptoTrial));
    [IntanBehaviour.missTrace(1:firstOptoTrial).opto] = a{:};
    [IntanBehaviour.MIFATrace(1:firstOptoTrial).opto] = a{:};
end

%% Separating trials 

% Hits 
optoTrials = find(cell2mat(arrayfun(@(s) s.opto,IntanBehaviour.cueHitTrace,'UniformOutput',false)) == 1);

IntanBehaviourOpto.cueHitTrace = IntanBehaviour.cueHitTrace(optoTrials);
IntanBehaviourOpto.hitTrace = IntanBehaviour.hitTrace(optoTrials);
IntanBehaviourOpto.MIHitTrace = IntanBehaviour.MIHitTrace(optoTrials);
IntanBehaviourOpto.reactionTime = IntanBehaviour.reactionTime(optoTrials);

IntanBehaviour.cueHitTrace(optoTrials) = [];
IntanBehaviour.hitTrace(optoTrials) = [];
IntanBehaviour.MIHitTrace(optoTrials) = [];
IntanBehaviour.reactionTime(optoTrials) = [];

if wavesPassed == 1
    WavesOpto.wavesHit = Waves.wavesHit(optoTrials);
    WavesOpto.wavesHitReward = Waves.wavesHitReward(optoTrials);
    WavesOpto.wavesMIHit = Waves.wavesMIHit(optoTrials);
    Waves.wavesHit(optoTrials) = [];
    Waves.wavesHitReward(optoTrials) = [];
    Waves.wavesMIHit(optoTrials) = [];
end

% Miss 
optoTrials = find(cell2mat(arrayfun(@(s) s.opto,IntanBehaviour.cueMissTrace,'UniformOutput',false)) == 1);
IntanBehaviourOpto.cueMissTrace = IntanBehaviour.cueMissTrace(optoTrials);
IntanBehaviour.cueMissTrace(optoTrials) = [];

if wavesPassed == 1
    WavesOpto.wavesMiss = Waves.wavesMiss(optoTrials);
    Waves.wavesMiss(optoTrials) = [];
end

% FA and MIFA
optoTrials = find(cell2mat(arrayfun(@(s) s.opto,IntanBehaviour.MIFATrace,'UniformOutput',false)) == 1);
IntanBehaviourOpto.MIFATrace = IntanBehaviour.MIFATrace(optoTrials);
IntanBehaviourOpto.missTrace = IntanBehaviour.missTrace(optoTrials);

IntanBehaviour.MIFATrace(optoTrials) = [];
IntanBehaviour.missTrace(optoTrials) = [];

if wavesPassed == 1
    WavesOpto.wavesFA = Waves.wavesFA(optoTrials);
    WavesOpto.wavesMIFA = Waves.wavesMIFA(optoTrials);
    Waves.wavesFA(optoTrials) = [];
    Waves.wavesMIFA(optoTrials) = [];
end

end

