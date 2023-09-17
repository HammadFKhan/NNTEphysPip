function [Behaviour] = readLever(parameters,lfpTime)

% if ~exist('parameters.experiment','var')
%     parameters.experiment = 'self';
%     disp('No experiment argument passed. Experiment type set to self initiated');
% end

if strcmp(parameters.experiment,'cue')
    cue = 1;
    disp('Experiment type set to cue initiated. . . ')
else
    cue = 0;
    disp('Experiment type set to self initiated. . . ')
end

if exist("lfpTime",'var')
    expFlag = 1;
    disp('Intan time data passed. Function set to experiment.');
else
    expFlag = 0;
    disp('Intan time data not passed. Function set to training.');
end
%% Reading file from arduino 
[enfile,enpath] = uigetfile('*.csv');
if isequal(enfile,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(enpath,enfile)]);
end

resting_position = 241;
flip = 1;
nlengthBeforePull = round(parameters.windowBeforePull/parameters.ts);
nlength = round(parameters.windowBeforePull/parameters.ts + parameters.windowAfterPull/parameters.ts + 1);
nlengthBeforeCue = round(parameters.windowBeforeCue/parameters.ts);
nlengthCue = round(parameters.windowBeforeCue/parameters.ts + parameters.windowAfterCue/parameters.ts + 1);

B = readmatrix([enpath,'/',enfile]);
Behaviour.leverTrace = (B(2:end,1) - resting_position)*flip;
Behaviour.time = (B(2:end,2) - B(2,2))/1e6; % time in seconds
Behaviour.nHit = B(end,3);
Behaviour.nMiss = B(end,4);
if cue==1 
    Behaviour.nCue = B(end,5); 
    Behaviour.nCueHit = Behaviour.nHit;
    Behaviour.nCueMiss = Behaviour.nCue-Behaviour.nCueHit;
end   
Behaviour.B = B(2:end,:);

%% Getting hit and miss timings
hitIndex = find(diff(B(:,3)) == 1);
hitTime = Behaviour.time(hitIndex);
if expFlag == 1
    hitLFPIndex = zeros(Behaviour.nHit,1);
    hitLFPTime = zeros(Behaviour.nHit,1);
    for i=1:Behaviour.nHit
        hitLFPIndex(i) = max(find(lfpTime<=hitTime(i)));
        hitLFPTime(i) = lfpTime(hitLFPIndex(i));
    end
    Behaviour.hit = [hitIndex hitTime hitLFPIndex hitLFPTime];
else
    Behaviour.hit = [hitIndex hitTime];
end



missIndex = find(diff(B(:,4)) == 1);
missTime = Behaviour.time(missIndex);
if expFlag == 1
    missLFPIndex = zeros(Behaviour.nMiss,1);
    missLFPTime = zeros(Behaviour.nMiss,1);
    for i=1:Behaviour.nMiss
        missLFPIndex(i) = max(find(lfpTime<=missTime(i)));
        missLFPTime(i) = lfpTime(missLFPIndex(i));
    end
    Behaviour.miss = [missIndex missTime missLFPIndex missLFPTime];
else
    Behaviour.miss = [missIndex missTime];
end
%% Getting cues, hit cues and miss cues
if cue == 1
    cueIndex = find(Behaviour.B(:,end) == 1);
    cueTime = Behaviour.time(cueIndex);
    if expFlag == 1
        cueLFPIndex = zeros(Behaviour.nCue,1);
        cueLFPTime = zeros(Behaviour.nCue,1);
        for i=1:Behaviour.nCue
            cueLFPIndex(i) = max(find(lfpTime<=cueTime(i)));
            cueLFPTime(i) = lfpTime(cueLFPIndex(i));
        end
        Behaviour.cue = [cueIndex cueTime cueLFPIndex cueLFPTime];
    else
        Behaviour.cue = [cueIndex cueTime];
    end
    
    % Getting cueHits and cueMisses 
    cueHitIndex = zeros(Behaviour.nHit,1);
    cueHitTime = zeros(Behaviour.nHit,1);
    cueHitPullIndex = zeros(Behaviour.nHit,1);
    cueHitPullTime = zeros(Behaviour.nHit,1);
    if expFlag == 1
        cueHitLFPIndex = zeros(Behaviour.nHit,1);
        cueHitLFPTime = zeros(Behaviour.nHit,1);
        cueHitPullLFPIndex = zeros(Behaviour.nHit,1);
        cueHitPullLFPTime = zeros(Behaviour.nHit,1);
    end

    a = zeros(Behaviour.nHit,1);
    
    for i=1:Behaviour.nHit
        a(i) = max(find(cueTime<hitTime(i))); % No need to check for reaction time. If there is a hit, there is a cue
        cueHitIndex(i) = cueIndex(a(i));
        cueHitTime(i) = cueTime(a(i));
        cueHitPullIndex(i) = hitIndex(i);
        cueHitPullTime(i) = hitTime(i);
        if expFlag == 1
            cueHitLFPIndex(i) = cueLFPIndex(a(i));
            cueHitLFPTime(i) = cueLFPTime(a(i));
            cueHitPullLFPIndex(i) = hitLFPIndex(i);
            cueHitPullLFPTime(i) = hitLFPTime(i);
        end
    end
    
    if expFlag == 1
        Behaviour.cueHit = [cueHitIndex cueHitTime cueHitLFPIndex cueHitLFPTime cueHitPullIndex cueHitPullTime cueHitPullLFPIndex cueHitPullLFPTime];
        Behaviour.reactionTime = cueHitPullLFPTime - cueHitLFPTime;
        Behaviour.meanReactionTime = mean(Behaviour.reactionTime,'all');
    else
        Behaviour.cueHit = [cueHitIndex cueHitTime cueHitPullIndex cueHitPullTime];
        Behaviour.reactionTime = cueHitPullTime - cueHitTime;
        Behaviour.meanReactionTime = mean(Behaviour.reactionTime,'all');
    end
    
    % For cue Miss trials
    Behaviour.cueMiss = Behaviour.cue;
    Behaviour.cueMiss(a,:) = [];
end

%% get lever traces for hits and miss 
st_hit1 = max(find(Behaviour.time < Behaviour.hit(1,2)-parameters.windowBeforePull));
if isempty(st_hit1)
    disp('First hit rejected');
    Behaviour.nHit = Behaviour.nHit-1;
    Behaviour.hit(1,:) = [];
end 
sp_hitend =  max(find(Behaviour.time < Behaviour.hit(end,2)+parameters.windowAfterPull));
if isempty(sp_hitend) 
    disp('Last hit rejected')
    Behaviour.nHit = Behaviour.nHit-1;
    Behaviour.hit(end,:) = [];
end

for i=1:Behaviour.nHit
    Behaviour.hitTrace(i).i1 = max(find(Behaviour.time < Behaviour.hit(i,2)-parameters.windowBeforePull));
    Behaviour.hitTrace(i).i0 = Behaviour.hit(i,1);
    Behaviour.hitTrace(i).i2 = max(find(Behaviour.time < Behaviour.hit(i,2)+parameters.windowAfterPull));
    Behaviour.hitTrace(i).rawtrace = Behaviour.leverTrace(Behaviour.hitTrace(i).i1:Behaviour.hitTrace(i).i2);
    Behaviour.hitTrace(i).rawtime = Behaviour.time(Behaviour.hitTrace(i).i1:Behaviour.hitTrace(i).i2) - Behaviour.time(Behaviour.hitTrace(i).i1);
    Behaviour.hitTrace(i).time1 = Behaviour.time(Behaviour.hitTrace(i).i1:Behaviour.hitTrace(i).i2);
    Behaviour.hitTrace(i).t1 = Behaviour.time(Behaviour.hitTrace(i).i1);
    Behaviour.hitTrace(i).t0 = Behaviour.hit(i,2);
    Behaviour.hitTrace(i).t2 = Behaviour.time(Behaviour.hitTrace(i).i2);
    if expFlag == 1
        [Behaviour.hitTrace(i).trace,Behaviour.hitTrace(i).time] = resample(Behaviour.hitTrace(i).rawtrace,Behaviour.hitTrace(i).rawtime,parameters.Fs,'spline');
        if (size(Behaviour.hitTrace(i).trace,1)<nlength)
            n1 = size(Behaviour.hitTrace(i).trace,1);
            n2 = nlength;
            Behaviour.hitTrace(i).trace = [Behaviour.hitTrace(i).trace;interp1(1:1:n1,Behaviour.hitTrace(i).trace,[n1+1:1:n2],'linear','extrap')'];
            Behaviour.hitTrace(i).time = [Behaviour.hitTrace(i).time;interp1(1:1:n1,Behaviour.hitTrace(i).time,[n1+1:1:n2],'linear','extrap')'];
        elseif (size(Behaviour.hitTrace(i).trace,1)>nlength)
            Behaviour.hitTrace(i).trace(nlength+1:end) = [];
            Behaviour.hitTrace(i).time(nlength+1:end) = []; 
        end
        Behaviour.hitTrace(i).LFPTime = Behaviour.hitTrace(i).time + (Behaviour.hit(i,4)-(nlengthBeforePull*parameters.ts));
        Behaviour.hitTrace(i).LFPIndex = ([Behaviour.hit(i,3)-nlengthBeforePull:1:nlengthBeforePull+Behaviour.hit(i,3)])';
    end
end

st_miss1 = max(find(Behaviour.time < Behaviour.miss(1,2)-parameters.windowBeforePull));
if isempty(st_miss1)
    disp('First miss rejected');
    Behaviour.nMiss = Behaviour.nMiss-1;
    Behaviour.miss(1,:) = [];
end 
sp_missend =  max(find(Behaviour.time < Behaviour.miss(end,2)+parameters.windowAfterPull));
if isempty(sp_missend)
    disp('Last miss rejected')
    Behaviour.nMiss = Behaviour.nMiss-1;
    Behaviour.miss(end,:) = [];
end

for i=1:Behaviour.nMiss
    Behaviour.missTrace(i).i1 = max(find(Behaviour.time < Behaviour.miss(i,2)-parameters.windowBeforePull));
    Behaviour.missTrace(i).i0 = Behaviour.miss(i,1);
    Behaviour.missTrace(i).i2 = max(find(Behaviour.time < Behaviour.miss(i,2)+parameters.windowAfterPull));
    Behaviour.missTrace(i).rawtrace = Behaviour.leverTrace(Behaviour.missTrace(i).i1:Behaviour.missTrace(i).i2);
    Behaviour.missTrace(i).rawtime = Behaviour.time(Behaviour.missTrace(i).i1:Behaviour.missTrace(i).i2) - Behaviour.time(Behaviour.missTrace(i).i1);
    Behaviour.missTrace(i).time1 = Behaviour.time(Behaviour.missTrace(i).i1:Behaviour.missTrace(i).i2);
    Behaviour.missTrace(i).t1 = Behaviour.time(Behaviour.missTrace(i).i1);
    Behaviour.missTrace(i).t0 = Behaviour.miss(i,2);
    Behaviour.missTrace(i).t2 = Behaviour.time(Behaviour.missTrace(i).i2);
    if expFlag == 1
        [Behaviour.missTrace(i).trace,Behaviour.missTrace(i).time] = resample(Behaviour.missTrace(i).rawtrace,Behaviour.missTrace(i).rawtime,parameters.Fs,'spline');
        if (size(Behaviour.missTrace(i).trace,1)<nlength)
            n1 = size(Behaviour.missTrace(i).trace,1);
            n2 = nlength;
            Behaviour.missTrace(i).trace = [Behaviour.missTrace(i).trace;interp1(1:1:n1,Behaviour.missTrace(i).trace,[n1+1:1:n2],'linear','extrap')'];
            Behaviour.missTrace(i).time = [Behaviour.missTrace(i).time;interp1(1:1:n1,Behaviour.missTrace(i).time,[n1+1:1:n2],'linear','extrap')'];
        elseif (size(Behaviour.missTrace(i).trace,1)>nlength)
            Behaviour.missTrace(i).trace(nlength+1:end) = [];
            Behaviour.missTrace(i).time(nlength+1:end) = []; 
        end
        Behaviour.missTrace(i).LFPTime = Behaviour.missTrace(i).time + (Behaviour.miss(i,4)-(nlengthBeforePull*parameters.ts));
        Behaviour.missTrace(i).LFPIndex = ([Behaviour.miss(i,3)-nlengthBeforePull:1:nlengthBeforePull+Behaviour.miss(i,3)])';
    end
end

if cue == 1 
    st_cuehit1 = max(find(Behaviour.time < Behaviour.cueHit(1,2)-parameters.windowBeforeCue));
    if isempty(st_cuehit1)
        disp('First cue hit rejected');
        Behaviour.nCueHit = Behaviour.nCueHit-1;
        Behaviour.cueHit(1,:) = [];
    end 
    sp_cuehitend =  max(find(Behaviour.time < Behaviour.cueHit(end,2)+parameters.windowAfterCue));
    if isempty(sp_cuehitend)
        disp('Last cue hit rejected')
        Behaviour.nCueHit = Behaviour.nCueHit-1;
        Behaviour.cueHit(end,:) = [];
    end
    for i=1:Behaviour.nCueHit
        Behaviour.cueHitTrace(i).i1 = max(find(Behaviour.time < Behaviour.cueHit(i,2)-parameters.windowBeforeCue));
        Behaviour.cueHitTrace(i).i0 = Behaviour.cueHit(i,1);
        Behaviour.cueHitTrace(i).i2 = max(find(Behaviour.time < Behaviour.cueHit(i,2)+parameters.windowAfterCue));
        Behaviour.cueHitTrace(i).rawtrace = Behaviour.leverTrace(Behaviour.cueHitTrace(i).i1:Behaviour.cueHitTrace(i).i2);
        Behaviour.cueHitTrace(i).rawtime = Behaviour.time(Behaviour.cueHitTrace(i).i1:Behaviour.cueHitTrace(i).i2) - Behaviour.time(Behaviour.cueHitTrace(i).i1);
        Behaviour.cueHitTrace(i).time1 = Behaviour.time(Behaviour.cueHitTrace(i).i1:Behaviour.cueHitTrace(i).i2);
        Behaviour.cueHitTrace(i).t1 = Behaviour.time(Behaviour.cueHitTrace(i).i1);
        Behaviour.cueHitTrace(i).t0 = Behaviour.cueHit(i,2);
        Behaviour.cueHitTrace(i).t2 = Behaviour.time(Behaviour.cueHitTrace(i).i2);
        if expFlag == 1
            [Behaviour.cueHitTrace(i).trace,Behaviour.cueHitTrace(i).time] = resample(Behaviour.cueHitTrace(i).rawtrace,Behaviour.cueHitTrace(i).rawtime,parameters.Fs,'spline');
            if (size(Behaviour.cueHitTrace(i).trace,1)<nlengthCue)
                n1 = size(Behaviour.cueHitTrace(i).trace,1);
                n2 = nlengthCue;
                Behaviour.cueHitTrace(i).trace = [Behaviour.cueHitTrace(i).trace;interp1(1:1:n1,Behaviour.cueHitTrace(i).trace,[n1+1:1:n2],'linear','extrap')'];
                Behaviour.cueHitTrace(i).time = [Behaviour.cueHitTrace(i).time;interp1(1:1:n1,Behaviour.cueHitTrace(i).time,[n1+1:1:n2],'linear','extrap')'];
            elseif (size(Behaviour.cueHitTrace(i).trace,1)>nlength)
                Behaviour.cueHitTrace(i).trace(nlength+1:end) = [];
                Behaviour.cueHitTrace(i).time(nlength+1:end) = []; 
            end
            Behaviour.cueHitTrace(i).LFPTime = Behaviour.cueHitTrace(i).time + (Behaviour.cueHit(i,4)-(nlengthBeforeCue*parameters.ts));
            Behaviour.cueHitTrace(i).LFPIndex = ([Behaviour.cueHit(i,3)-nlengthBeforeCue:1:nlengthBeforeCue+Behaviour.cueHit(i,3)])';
        end
    end
    
    % Getting cueMiss traces
    if Behaviour.nCueMiss > 0
        st_cuemiss1 = max(find(Behaviour.time < Behaviour.cueMiss(1,2)-parameters.windowBeforeCue));
        if isempty(st_cuemiss1)
            disp('First cue miss rejected');
            Behaviour.nCueMiss = Behaviour.nCueMiss-1;
            Behaviour.cueMiss(1,:) = [];
        end 
        sp_cuemissend =  max(find(Behaviour.time < Behaviour.cueMiss(end,2)+parameters.windowAfterCue));
        if isempty(sp_cuemissend)
            disp('Last cue hit rejected')
            Behaviour.nCueMiss = Behaviour.nCueMiss-1;
            Behaviour.cueMiss(end,:) = [];
        end
        for i=1:Behaviour.nCueMiss
            Behaviour.cueMissTrace(i).i1 = max(find(Behaviour.time < Behaviour.cueMiss(i,2)-parameters.windowBeforeCue));
            Behaviour.cueMissTrace(i).i0 = Behaviour.cueMiss(i,1);
            Behaviour.cueMissTrace(i).i2 = max(find(Behaviour.time < Behaviour.cueMiss(i,2)+parameters.windowAfterCue));
            Behaviour.cueMissTrace(i).rawtrace = Behaviour.leverTrace(Behaviour.cueMissTrace(i).i1:Behaviour.cueMissTrace(i).i2);
            Behaviour.cueMissTrace(i).rawtime = Behaviour.time(Behaviour.cueMissTrace(i).i1:Behaviour.cueMissTrace(i).i2) - Behaviour.time(Behaviour.cueMissTrace(i).i1);
            Behaviour.cueMissTrace(i).time1 = Behaviour.time(Behaviour.cueMissTrace(i).i1:Behaviour.cueMissTrace(i).i2);
            Behaviour.cueMissTrace(i).t1 = Behaviour.time(Behaviour.cueMissTrace(i).i1);
            Behaviour.cueMissTrace(i).t0 = Behaviour.cueMiss(i,2);
            Behaviour.cueMissTrace(i).t2 = Behaviour.time(Behaviour.cueMissTrace(i).i2);
            if expFlag == 1
                [Behaviour.cueMissTrace(i).trace,Behaviour.cueMissTrace(i).time] = resample(Behaviour.cueMissTrace(i).rawtrace,Behaviour.cueMissTrace(i).rawtime,parameters.Fs,'spline');
                if (size(Behaviour.cueMissTrace(i).trace,1)<nlengthCue)
                    n1 = size(Behaviour.cueMissTrace(i).trace,1);
                    n2 = nlengthCue;
                    Behaviour.cueMissTrace(i).trace = [Behaviour.cueMissTrace(i).trace;interp1(1:1:n1,Behaviour.cueMissTrace(i).trace,[n1+1:1:n2],'linear','extrap')'];
                    Behaviour.cueMissTrace(i).time = [Behaviour.cueMissTrace(i).time;interp1(1:1:n1,Behaviour.cueMissTrace(i).time,[n1+1:1:n2],'linear','extrap')'];
                elseif (size(Behaviour.cueMissTrace(i).trace,1)>nlength)
                    Behaviour.cueMissTrace(i).trace(nlength+1:end) = [];
                    Behaviour.cueMissTrace(i).time(nlength+1:end) = []; 
                end
                Behaviour.cueMissTrace(i).LFPTime = Behaviour.cueMissTrace(i).time + (Behaviour.cueMiss(i,4)-(nlengthBeforeCue*parameters.ts));
                Behaviour.cueMissTrace(i).LFPIndex = ([Behaviour.cueMiss(i,3)-nlengthBeforeCue:1:nlengthBeforeCue+Behaviour.cueMiss(i,3)])';
            end
        end
    end
end
