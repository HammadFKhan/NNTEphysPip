function [IntanBehaviour] = readLeverIntan(parameters,lfpTime,leverTrace,rewardTrace)


intanFs = 10000;

resting_position = 242*5/1024;
flip = 1;
nlengthBeforePull = round(parameters.windowAfterPull/parameters.ts);
nlength = round(parameters.windowAfterPull/parameters.ts + parameters.windowAfterPull/parameters.ts + 1);

IntanBehaviour.leverTrace = resample(((double(leverTrace)- resting_position)*flip),parameters.Fs,intanFs);
IntanBehaviour.time = lfpTime; % time in seconds
IntanBehaviour.rewardTrace = downsample(rewardTrace,round(intanFs/parameters.Fs),2); 

rewardIndex = find(diff(IntanBehaviour.rewardTrace)==1)+1;

IntanBehaviour.nHit = size(rewardIndex,2);
% IntanBehaviour.nMiss = B(end,4);

%% Getting hit traces and timings

for i=1:IntanBehaviour.nHit
%     IntanBehaviour.hit(i) = [rewardIndex(i) lfpTime(rewardIndex(i)) rewardIndex(i) lfpTime(rewardIndex(i))];
    IntanBehaviour.hitTrace(i).trace = IntanBehaviour.leverTrace(rewardIndex(i)-parameters.windowBeforePull*parameters.Fs:rewardIndex(i)+parameters.windowAfterPull*parameters.Fs)';
    IntanBehaviour.hitTrace(i).time = IntanBehaviour.time(rewardIndex(i)-parameters.windowBeforePull*parameters.Fs:rewardIndex(i)+parameters.windowAfterPull*parameters.Fs)';
    IntanBehaviour.hitTrace(i).LFPIndex = ([rewardIndex(i)-parameters.windowBeforePull*parameters.Fs:1:rewardIndex(i)+parameters.windowAfterPull*parameters.Fs])';
    IntanBehaviour.hitTrace(i).LFPtime = IntanBehaviour.hitTrace(i).time;
end


% %% Gettting Hit and miss timings
% hitIndex = find(diff(B(:,3)) == 1) + 1;
% hitTime = IntanBehaviour.time(hitIndex);
% hitLFPIndex = zeros(IntanBehaviour.nHit,1);
% hitLFPTime = zeros(IntanBehaviour.nHit,1);
% for i=1:IntanBehaviour.nHit
%     hitLFPIndex(i) = max(find(lfpTime<hitTime(i)));
%     hitLFPTime(i) = lfpTime(hitLFPIndex(i));
% end
% IntanBehaviour.hit = [hitIndex hitTime hitLFPIndex hitLFPTime];
% 
% missIndex = find(diff(B(:,4)) == 1) + 1;
% missTime = IntanBehaviour.time(missIndex);
% missLFPIndex = zeros(IntanBehaviour.nMiss,1);
% missLFPTime = zeros(IntanBehaviour.nMiss,1);
% for i=1:IntanBehaviour.nMiss
%     missLFPIndex(i) = max(find(lfpTime<missTime(i)));
%     missLFPTime(i) = lfpTime(missLFPIndex(i));
% end
% IntanBehaviour.miss = [missIndex missTime missLFPIndex missLFPTime];
% 
% %% get lever traces for hits and miss 
% 
% for i=1:IntanBehaviour.nHit
%     IntanBehaviour.hitTrace(i).i1 = max(find(IntanBehaviour.time < IntanBehaviour.hit(i,2)-parameters.windowBeforePull));
%     IntanBehaviour.hitTrace(i).i0 = IntanBehaviour.hit(i,1);
%     IntanBehaviour.hitTrace(i).i2 = max(find(IntanBehaviour.time < IntanBehaviour.hit(i,2)+parameters.windowBeforePull));
%     IntanBehaviour.hitTrace(i).rawtrace = IntanBehaviour.leverTrace(IntanBehaviour.hitTrace(i).i1:IntanBehaviour.hitTrace(i).i2);
%     IntanBehaviour.hitTrace(i).rawtime = IntanBehaviour.time(IntanBehaviour.hitTrace(i).i1:IntanBehaviour.hitTrace(i).i2) - IntanBehaviour.time(IntanBehaviour.hitTrace(i).i1);
%     IntanBehaviour.hitTrace(i).time1 = IntanBehaviour.time(IntanBehaviour.hitTrace(i).i1:IntanBehaviour.hitTrace(i).i2);
%     IntanBehaviour.hitTrace(i).t1 = IntanBehaviour.time(IntanBehaviour.hitTrace(i).i1);
%     IntanBehaviour.hitTrace(i).t0 = IntanBehaviour.hit(i,2);
%     IntanBehaviour.hitTrace(i).t2 = IntanBehaviour.time(IntanBehaviour.hitTrace(i).i2);
%     [IntanBehaviour.hitTrace(i).trace,IntanBehaviour.hitTrace(i).time] = resample(IntanBehaviour.hitTrace(i).rawtrace,IntanBehaviour.hitTrace(i).rawtime,parameters.Fs,'spline');
%     if (size(IntanBehaviour.hitTrace(i).trace,1)<nlength)
%         n1 = size(IntanBehaviour.hitTrace(i).trace,1);
%         n2 = nlength;
%         IntanBehaviour.hitTrace(i).trace = [IntanBehaviour.hitTrace(i).trace;interp1(1:1:n1,IntanBehaviour.hitTrace(i).trace,[n1+1:1:n2],'linear','extrap')'];
%         IntanBehaviour.hitTrace(i).time = [IntanBehaviour.hitTrace(i).time;interp1(1:1:n1,IntanBehaviour.hitTrace(i).time,[n1+1:1:n2],'linear','extrap')'];
%     elseif (size(IntanBehaviour.hitTrace(i).trace,1)>nlength)
%         IntanBehaviour.hitTrace(i).trace(nlength+1:end) = [];
%         IntanBehaviour.hitTrace(i).time(nlength+1:end) = []; 
%     end
%     IntanBehaviour.hitTrace(i).LFPTime = IntanBehaviour.hitTrace(i).time + (IntanBehaviour.hit(i,4)-(nlengthBeforePull*parameters.ts));
%     IntanBehaviour.hitTrace(i).LFPIndex = ([IntanBehaviour.hit(i,3)-nlengthBeforePull:1:nlengthBeforePull+IntanBehaviour.hit(i,3)])';
% end
% 
% for i=1:IntanBehaviour.nMiss
%     IntanBehaviour.missTrace(i).i1 = max(find(IntanBehaviour.time < IntanBehaviour.miss(i,2)-parameters.windowBeforePull));
%     IntanBehaviour.missTrace(i).i0 = IntanBehaviour.miss(i,1);
%     IntanBehaviour.missTrace(i).i2 = max(find(IntanBehaviour.time < IntanBehaviour.miss(i,2)+parameters.windowBeforePull));
%     IntanBehaviour.missTrace(i).rawtrace = IntanBehaviour.leverTrace(IntanBehaviour.missTrace(i).i1:IntanBehaviour.missTrace(i).i2);
%     IntanBehaviour.missTrace(i).rawtime = IntanBehaviour.time(IntanBehaviour.missTrace(i).i1:IntanBehaviour.missTrace(i).i2) - IntanBehaviour.time(IntanBehaviour.missTrace(i).i1);
%     IntanBehaviour.missTrace(i).time1 = IntanBehaviour.time(IntanBehaviour.missTrace(i).i1:IntanBehaviour.missTrace(i).i2);
%     IntanBehaviour.missTrace(i).t1 = IntanBehaviour.time(IntanBehaviour.missTrace(i).i1);
%     IntanBehaviour.missTrace(i).t0 = IntanBehaviour.miss(i,2);
%     IntanBehaviour.missTrace(i).t2 = IntanBehaviour.time(IntanBehaviour.missTrace(i).i2);
%     [IntanBehaviour.missTrace(i).trace,IntanBehaviour.missTrace(i).time] = resample(IntanBehaviour.missTrace(i).rawtrace,IntanBehaviour.missTrace(i).rawtime,parameters.Fs,'spline');
%     if (size(IntanBehaviour.missTrace(i).trace,1)<nlength)
%         n1 = size(IntanBehaviour.missTrace(i).trace,1);
%         n2 = nlength;
%         IntanBehaviour.missTrace(i).trace = [IntanBehaviour.missTrace(i).trace;interp1(1:1:n1,IntanBehaviour.missTrace(i).trace,[n1+1:1:n2],'linear','extrap')'];
%         IntanBehaviour.missTrace(i).time = [IntanBehaviour.missTrace(i).time;interp1(1:1:n1,IntanBehaviour.missTrace(i).time,[n1+1:1:n2],'linear','extrap')'];
%     elseif (size(IntanBehaviour.missTrace(i).trace,1)>nlength)
%         IntanBehaviour.missTrace(i).trace(nlength+1:end) = [];
%         IntanBehaviour.missTrace(i).time(nlength+1:end) = []; 
%     end
%     IntanBehaviour.missTrace(i).LFPTime = IntanBehaviour.missTrace(i).time + (IntanBehaviour.miss(i,4)-(nlengthBeforePull*parameters.ts));
%     IntanBehaviour.missTrace(i).LFPIndex = ([IntanBehaviour.miss(i,3)-nlengthBeforePull:1:nlengthBeforePull+IntanBehaviour.miss(i,3)])';
% end



