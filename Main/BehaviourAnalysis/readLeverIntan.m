function [IntanBehaviour] = readLeverIntan(parameters,lfpTime,leverTrace,dig_in_data,Behaviour)

if strcmp(parameters.experiment,'self')
    rewardTrace = dig_in_data(1,:);
    disp('Experiment set to self initiated.');
    cue = 0;
elseif strcmp(parameters.experiment,'cue')
    rewardTrace = dig_in_data(1,:);
    cueTrace = dig_in_data(2,:);
    disp('Experiment set to cue initiated.');
    cue = 1;
end

if parameters.opto == 1
    optoTrace = dig_in_data(3,:);
    disp('Experiment has opto trials.');
else
    disp('Experiment has no opto trials.');
end

intanFs = 5000;

resting_position = 241*5/1024;
flip = 1;
nlengthBeforePull = round(parameters.windowAfterPull/parameters.ts);
nlength = round(parameters.windowAfterPull/parameters.ts + parameters.windowAfterPull/parameters.ts + 1);
nlengthBeforeCue = round(parameters.windowBeforeCue/parameters.ts);
nlengthCue = round(parameters.windowBeforeCue/parameters.ts + parameters.windowAfterCue/parameters.ts + 1);

IntanBehaviour.leverTrace = resample(((double(leverTrace)- resting_position)*flip),parameters.Fs,intanFs);
IntanBehaviour.time = lfpTime; % time in seconds
IntanBehaviour.rewardTrace = downsample(rewardTrace,round(intanFs/parameters.Fs),2);
rewardIndex = find(diff(IntanBehaviour.rewardTrace)==1)+1;
if cue == 1
    IntanBehaviour.cueTrace = downsample(cueTrace,round(intanFs/parameters.Fs),2); 
    cueIndex = find(diff(IntanBehaviour.cueTrace)==1)+1;
end
if parameters.opto == 1
    IntanBehaviour.optoTrace = downsample(optoTrace,round(intanFs/parameters.Fs),2);
    optoIndex = find(diff(IntanBehaviour.optoTrace)==1)+1;
end

IntanBehaviour.nCueHit = size(rewardIndex,2);
IntanBehaviour.nCueMiss = Behaviour.nCueMiss;
IntanBehaviour.nOpto = size(optoIndex,2);

% Estimating the threshold for reward
IntanBehaviour.threshold = mean(IntanBehaviour.leverTrace(rewardIndex),'all');

%% Getting hit traces and timings
st_hit1 = rewardIndex(1)-parameters.windowBeforePull*parameters.Fs;
if st_hit1 <= 0
    disp('First cue hit rejected');
    IntanBehaviour.nCueHit = Behaviour.nCueHit-1;
    rewardIndex(1) = [];
end
sp_hitend = rewardIndex(end)+parameters.windowAfterPull*parameters.Fs;
if  sp_hitend > length(IntanBehaviour.time)
    disp('Last cue hit rejected');
    IntanBehaviour.nCueHit = Behaviour.nCueHit-1;
    rewardIndex(end) = [];
end 

for i=1:IntanBehaviour.nCueHit
%     IntanBehaviour.hit(i) = [rewardIndex(i) lfpTime(rewardIndex(i)) rewardIndex(i) lfpTime(rewardIndex(i))];
    IntanBehaviour.hitTrace(i).trace = IntanBehaviour.leverTrace(rewardIndex(i)-parameters.windowBeforePull*parameters.Fs:rewardIndex(i)+parameters.windowAfterPull*parameters.Fs)';
    IntanBehaviour.hitTrace(i).time = (0:1/parameters.Fs:(size(IntanBehaviour.hitTrace(i).trace,1)-1)*1/parameters.Fs)';
    IntanBehaviour.hitTrace(i).LFPIndex = ([rewardIndex(i)-parameters.windowBeforePull*parameters.Fs:1:rewardIndex(i)+parameters.windowAfterPull*parameters.Fs])';
    IntanBehaviour.hitTrace(i).LFPtime = IntanBehaviour.time(rewardIndex(i)-parameters.windowBeforePull*parameters.Fs:rewardIndex(i)+parameters.windowAfterPull*parameters.Fs)';
end


%% Getting miss traces and timings

% figure();plot(1:1:10001,(5/1024)*Behaviour.leverTrace(Behaviour.miss(15,1)-5000:Behaviour.miss(15,1)+5000));hold on;
% plot(0.1:0.1:10000.1,IntanBehaviour.leverTrace(Behaviour.miss(15,3)-50000:Behaviour.miss(15,3)+50000));
correctionWindow = 800; % in number of points in LFPFs
tol = 0.005;
nDiffSlope = 10;
disp('Finding miss trials in the Intan data ...');
for i=1:Behaviour.nMiss
    missIndexAr = Behaviour.miss(i,3);
    trace1 = IntanBehaviour.leverTrace(missIndexAr-correctionWindow:missIndexAr+correctionWindow);
    misstrigs1 = find(trace1 <IntanBehaviour.threshold+tol & trace1>IntanBehaviour.threshold-tol);  
    if isempty(misstrigs1)
        missIndex(i) = NaN;
        disp('Suss miss trial removed');
        continue;
    end
    % Checking slope 
    for j=1:size(misstrigs1,2)
        % Checking edge cases, rejects all the edge cases 
        if((misstrigs1(j)+nDiffSlope >= (correctionWindow*2+1)) || (misstrigs1(j)-nDiffSlope <= 0)) 
            misstrigs1(j) = NaN;
            continue
        end
        % Checking slope, reject all negative slope 
        slope = mean(trace1(misstrigs1(j):misstrigs1(j)+nDiffSlope)) - mean(trace1(misstrigs1(j)-nDiffSlope:misstrigs1(j)));
        if slope < 0
            misstrigs1(j) = NaN;
        end
    end
    misstrigs1 = misstrigs1(~isnan(misstrigs1));
    missIndex(i) =  missIndexAr - correctionWindow + min(misstrigs1);
end

missIndex = removeNaNRows(missIndex');

IntanBehaviour.nMiss = size(missIndex,1);

for i=1:IntanBehaviour.nMiss
%     IntanBehaviour.hit(i) = [rewardIndex(i) lfpTime(rewardIndex(i)) rewardIndex(i) lfpTime(rewardIndex(i))];
    IntanBehaviour.missTrace(i).trace = IntanBehaviour.leverTrace(missIndex(i)-parameters.windowBeforePull*parameters.Fs:missIndex(i)+parameters.windowAfterPull*parameters.Fs)';
    IntanBehaviour.missTrace(i).time = (0:1/parameters.Fs:(size(IntanBehaviour.hitTrace(i).trace,1)-1)*1/parameters.Fs)';
    IntanBehaviour.missTrace(i).LFPIndex = ([missIndex(i)-parameters.windowBeforePull*parameters.Fs:1:missIndex(i)+parameters.windowAfterPull*parameters.Fs])';
    IntanBehaviour.missTrace(i).LFPtime = IntanBehaviour.time(missIndex(i)-parameters.windowBeforePull*parameters.Fs:missIndex(i)+parameters.windowAfterPull*parameters.Fs)';
end

%% Getting cue Hit traces
a = zeros(IntanBehaviour.nCueHit,1);
shortTrialFlag = 0; % 0 - none, 1 - first, 2 = end , 3 = both

for i=1:IntanBehaviour.nCueHit
    a(i) = max(find(cueIndex<=rewardIndex(i)));
    IntanBehaviour.cueHit(i,1) = cueIndex(a(i)); % Cue index 
    IntanBehaviour.cueHit(i,2) = rewardIndex(i); % pull index for hits
    if ~isempty(find(optoIndex == IntanBehaviour.cueHit(i,1))) 
        IntanBehaviour.cueHitTrace(i).opto = 1;
    else 
        IntanBehaviour.cueHitTrace(i).opto = 0;
    end 
    st = IntanBehaviour.cueHit(i,1)-parameters.windowBeforeCue*parameters.Fs;
    sp = IntanBehaviour.cueHit(i,1)+parameters.windowAfterCue*parameters.Fs;
    if st <= 0
        st = 1;
        disp('Short Cue Hit trial Detected');
        shortTrialFlag = 1;
    end
    if sp >= size(IntanBehaviour.leverTrace,2)
        sp = size(IntanBehaviour.leverTrace,2);
        disp('Short Cue Hit trial Detected');
        if shortTrialFlag == 1 shortTrialFlag = 3, else shortTrialFlag = 2, end
    end 
    
    IntanBehaviour.cueHitTrace(i).trace = IntanBehaviour.leverTrace(st:sp)';
    IntanBehaviour.cueHitTrace(i).time = (0:1/parameters.Fs:(size(IntanBehaviour.cueHitTrace(i).trace,1)-1)*1/parameters.Fs)' - parameters.windowBeforeCue;
    IntanBehaviour.cueHitTrace(i).LFPIndex = ([st:1:sp])';
    IntanBehaviour.cueHitTrace(i).LFPtime = IntanBehaviour.time(st:sp)';
end

if shortTrialFlag == 1
    IntanBehaviour.cueHit(1,:) = [];
    IntanBehaviour.cueHitTrace(1) = [];
    IntanBehaviour.nCueHit = IntanBehaviour.nCueHit-1;
elseif shortTrialFlag == 2
    IntanBehaviour.cueHit(end,:) = [];
    IntanBehaviour.cueHitTrace(end) = [];
    IntanBehaviour.nCueHit = IntanBehaviour.nCueHit-1;
elseif shortTrialFlag == 3
    IntanBehaviour.cueHit(1,:) = [];
    IntanBehaviour.cueHit(end,:) = [];
    IntanBehaviour.cueHitTrace(1) = [];
    IntanBehaviour.cueHitTrace(end) = [];
    IntanBehaviour.nCueHit = IntanBehaviour.nCueHit-2;
end

% Reaction Time 
IntanBehaviour.reactionTime = diff(IntanBehaviour.cueHit,1,2)/parameters.Fs; % in seconds

%% Getting cue Miss traces
IntanBehaviour.cueMiss = cueIndex';
IntanBehaviour.cueMiss(a) = [];
IntanBehaviour.nCueMiss = size(IntanBehaviour.cueMiss,1);

st1 = IntanBehaviour.cueMiss(1,1)-parameters.windowBeforeCue*parameters.Fs;
spend = IntanBehaviour.cueMiss(end,1)+parameters.windowAfterCue*parameters.Fs;
if st1 <= 0
    disp('Short Cue Miss trial Detected');
    IntanBehaviour.cueMiss(1) = [];
    IntanBehaviour.nCueMiss = IntanBehaviour.nCueMiss-1;
end
if spend >= size(IntanBehaviour.leverTrace,2)
    disp('Short Cue Miss trial Detected');
    IntanBehaviour.cueMiss(end) = [];
    IntanBehaviour.nCueMiss = IntanBehaviour.nCueMiss-1;
end

for i=1:IntanBehaviour.nCueMiss
    if ~isempty(find(optoIndex == IntanBehaviour.cueMiss(i,1))) 
        IntanBehaviour.cueMissTrace(i).opto = 1;
    else
        IntanBehaviour.cueMissTrace(i).opto = 0;
    end 
    st = IntanBehaviour.cueMiss(i,1)-parameters.windowBeforeCue*parameters.Fs;
    sp = IntanBehaviour.cueMiss(i,1)+parameters.windowAfterCue*parameters.Fs;
    IntanBehaviour.cueMissTrace(i).trace = IntanBehaviour.leverTrace(st:sp)';
    IntanBehaviour.cueMissTrace(i).time = (0:1/parameters.Fs:(size(IntanBehaviour.cueMissTrace(i).trace,1)-1)*1/parameters.Fs)' - parameters.windowBeforeCue;
    IntanBehaviour.cueMissTrace(i).LFPIndex = ([st:1:sp])';
    IntanBehaviour.cueMissTrace(i).LFPtime = IntanBehaviour.time(st:sp)';
end