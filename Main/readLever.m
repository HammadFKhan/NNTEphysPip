function [Behaviour] = readLever(parameters,lfpTime)
%% Reading file from arduino 
[enfile,enpath] = uigetfile('*.csv');
if isequal(enfile,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(enpath,enfile)]);
end

resting_position = 245;
flip = 1;
nlengthBeforePull = round(parameters.windowAfterPull/parameters.ts);
nlength = round(parameters.windowAfterPull/parameters.ts + parameters.windowAfterPull/parameters.ts + 1);

B = readmatrix([enpath,'/',enfile]);
Behaviour.leverTrace = (B(2:end,1) - resting_position)*flip;
Behaviour.time = (B(2:end,2) - B(2,2))/1e6; % time in seconds
Behaviour.nHit = B(end,3);
Behaviour.nMiss = B(end,4);
Behaviour.B = B(2:end,:);

%% Gettting Hit and miss timings
hitIndex = find(diff(B(:,3)) == 1) + 1;
hitTime = Behaviour.time(hitIndex);
hitLFPIndex = zeros(Behaviour.nHit,1);
hitLFPTime = zeros(Behaviour.nHit,1);
for i=1:Behaviour.nHit
    hitLFPIndex(i) = max(find(lfpTime<hitTime(i)));
    hitLFPTime(i) = lfpTime(hitLFPIndex(i));
end
Behaviour.hit = [hitIndex hitTime hitLFPIndex hitLFPTime];

missIndex = find(diff(B(:,4)) == 1) + 1;
missTime = Behaviour.time(missIndex);
missLFPIndex = zeros(Behaviour.nMiss,1);
missLFPTime = zeros(Behaviour.nMiss,1);
for i=1:Behaviour.nMiss
    missLFPIndex(i) = max(find(lfpTime<missTime(i)));
    missLFPTime(i) = lfpTime(missLFPIndex(i));
end
Behaviour.miss = [missIndex missTime missLFPIndex missLFPTime];

%% get lever traces for hits and miss 

for i=1:Behaviour.nHit
    Behaviour.hitTrace(i).i1 = max(find(Behaviour.time < Behaviour.hit(i,2)-parameters.windowBeforePull));
    Behaviour.hitTrace(i).i0 = Behaviour.hit(i,1);
    Behaviour.hitTrace(i).i2 = max(find(Behaviour.time < Behaviour.hit(i,2)+parameters.windowBeforePull));
    Behaviour.hitTrace(i).rawtrace = Behaviour.leverTrace(Behaviour.hitTrace(i).i1:Behaviour.hitTrace(i).i2);
    Behaviour.hitTrace(i).rawtime = Behaviour.time(Behaviour.hitTrace(i).i1:Behaviour.hitTrace(i).i2) - Behaviour.time(Behaviour.hitTrace(i).i1);
    Behaviour.hitTrace(i).time1 = Behaviour.time(Behaviour.hitTrace(i).i1:Behaviour.hitTrace(i).i2);
    Behaviour.hitTrace(i).t1 = Behaviour.time(Behaviour.hitTrace(i).i1);
    Behaviour.hitTrace(i).t0 = Behaviour.hit(i,2);
    Behaviour.hitTrace(i).t2 = Behaviour.time(Behaviour.hitTrace(i).i2);
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

for i=1:Behaviour.nMiss
    Behaviour.missTrace(i).i1 = max(find(Behaviour.time < Behaviour.miss(i,2)-parameters.windowBeforePull));
    Behaviour.missTrace(i).i0 = Behaviour.miss(i,1);
    Behaviour.missTrace(i).i2 = max(find(Behaviour.time < Behaviour.miss(i,2)+parameters.windowBeforePull));
    Behaviour.missTrace(i).rawtrace = Behaviour.leverTrace(Behaviour.missTrace(i).i1:Behaviour.missTrace(i).i2);
    Behaviour.missTrace(i).rawtime = Behaviour.time(Behaviour.missTrace(i).i1:Behaviour.missTrace(i).i2) - Behaviour.time(Behaviour.missTrace(i).i1);
    Behaviour.missTrace(i).time1 = Behaviour.time(Behaviour.missTrace(i).i1:Behaviour.missTrace(i).i2);
    Behaviour.missTrace(i).t1 = Behaviour.time(Behaviour.missTrace(i).i1);
    Behaviour.missTrace(i).t0 = Behaviour.miss(i,2);
    Behaviour.missTrace(i).t2 = Behaviour.time(Behaviour.missTrace(i).i2);
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




