function [Behaviour] = readLever2(parameters)
%% Reading file from arduino 
[enfile,enpath] = uigetfile('*.csv');
if isequal(enfile,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(enpath,enfile)]);
end
resting_position = 242;
flip = 1;
nlengthBeforePull = round(parameters.windowAfterPull/parameters.ts);
nlength = round(parameters.windowAfterPull/parameters.ts + parameters.windowAfterPull/parameters.ts + 1);
​
B = readmatrix([enpath,'/',enfile]);
Behaviour.leverTrace = (B(2:end,1) - resting_position)*flip;
Behaviour.time = (B(2:end,2) - B(2,2))/1e6; % time in seconds
Behaviour.nHit = B(end,3);
Behaviour.nMiss = B(end,4);
Behaviour.B = B(2:end,:);
​
%% Gettting Hit and miss timings
hitIndex = find(diff(B(:,3)) == 1) + 1;
hitTime = Behaviour.time(hitIndex);
hitLFPIndex = zeros(Behaviour.nHit,1);
hitLFPTime = zeros(Behaviour.nHit,1);
Behaviour.hit = [hitIndex hitTime];
​
missIndex = find(diff(B(:,4)) == 1) + 1;
missTime = Behaviour.time(missIndex);
missLFPIndex = zeros(Behaviour.nMiss,1);
missLFPTime = zeros(Behaviour.nMiss,1);
Behaviour.miss = [missIndex missTime];
​
%% get lever traces for hits and miss 
​
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
end
​
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
end
​