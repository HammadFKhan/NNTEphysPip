function [IntanBehaviour] = readLeverIntan(parameters,lfpTime,leverTrace,dig_in_data,Behaviour,plotOption)

if strcmp(parameters.experiment,'self')
    rewardTrace = dig_in_data(2,:);
    disp('Experiment set to self initiated.');
    cue = 0;
elseif strcmp(parameters.experiment,'cue')
    if size(dig_in_data,1)<3
        rewardTrace = dig_in_data(1,:);
        cueTrace = dig_in_data(2,:);
    else
        rewardTrace = dig_in_data(2,:);
        cueTrace = dig_in_data(3,:);
    end
    disp('Experiment set to cue initiated.');
    cue = 1;
end

if parameters.opto == 1
    optoTrace = dig_in_data(4,:);
    disp('Experiment has opto trials.');
else
    disp('Experiment has no opto trials.');
end

intanFs = parameters.IntanFs;

resting_position = 285*5/1024;
flip = 1;
nlengthBeforePull = round(parameters.windowAfterPull/parameters.ts);
nlength = round(parameters.windowAfterPull/parameters.ts + parameters.windowAfterPull/parameters.ts + 1);
nlengthBeforeCue = round(parameters.windowBeforeCue/parameters.ts);
nlengthCue = round(parameters.windowBeforeCue/parameters.ts + parameters.windowAfterCue/parameters.ts + 1);

IntanBehaviour.leverTrace = resample(((double(leverTrace)- resting_position)*flip),parameters.Fs,intanFs);
IntanBehaviour.time = downsample(lfpTime,round(intanFs/parameters.Fs),1); % time in seconds
IntanBehaviour.rewardTrace = downsample(rewardTrace,round(intanFs/parameters.Fs),1);
rewardIndex = find(diff(IntanBehaviour.rewardTrace)==1)+1;
if cue == 1
    IntanBehaviour.cueTrace = downsample(cueTrace,round(intanFs/parameters.Fs),1); 
    cueIndex = find(diff(IntanBehaviour.cueTrace)==1)+1;
    IntanBehaviour.nCueHit = size(rewardIndex,2);
    IntanBehaviour.nCueMiss = Behaviour.nCueMiss;
end
if parameters.opto == 1
    IntanBehaviour.optoTrace = downsample(optoTrace,round(intanFs/parameters.Fs),1);
    optoIndex = find(diff(IntanBehaviour.optoTrace)==1)+1;
    IntanBehaviour.nOpto = size(optoIndex,2);
end

IntanBehaviour.nHit = size(rewardIndex,2);


% Estimating the threshold for reward
IntanBehaviour.threshold = mean(IntanBehaviour.leverTrace(rewardIndex),'all');

%% Getting hit traces and timings

st_hit1 = rewardIndex(1)-parameters.windowBeforePull*parameters.Fs;
if st_hit1 <= 0
    disp('First cue hit rejected');
    IntanBehaviour.nHit = IntanBehaviour.nHit-1;
    if cue == 1
        IntanBehaviour.nCueHit = IntanBehaviour.nCueHit-1;
    end
    rewardIndex(1) = [];
end
sp_hitend = rewardIndex(end)+parameters.windowAfterPull*parameters.Fs;
if  sp_hitend > length(IntanBehaviour.time)
    disp('Last cue hit rejected');
    IntanBehaviour.nHit = Behaviour.nHit-1;
    if cue == 1
        IntanBehaviour.nCueHit = IntanBehaviour.nCueHit-1;
    end
    rewardIndex(end) = [];
end 

for i=1:IntanBehaviour.nHit
%     IntanBehaviour.hit(i) = [rewardIndex(i) lfpTime(rewardIndex(i)) rewardIndex(i) lfpTime(rewardIndex(i))];
    IntanBehaviour.hitTrace(i).trace = IntanBehaviour.leverTrace(rewardIndex(i)-parameters.windowBeforePull*parameters.Fs:rewardIndex(i)+parameters.windowAfterPull*parameters.Fs)';
    IntanBehaviour.hitTrace(i).time = (0:1/parameters.Fs:(size(IntanBehaviour.hitTrace(i).trace,1)-1)*1/parameters.Fs)' - parameters.windowBeforePull;
    IntanBehaviour.hitTrace(i).LFPIndex = ([rewardIndex(i)-parameters.windowBeforePull*parameters.Fs:1:rewardIndex(i)+parameters.windowAfterPull*parameters.Fs])';
    IntanBehaviour.hitTrace(i).LFPtime = IntanBehaviour.time(rewardIndex(i)-parameters.windowBeforePull*parameters.Fs:rewardIndex(i)+parameters.windowAfterPull*parameters.Fs)';
end

%% Getting cue Hit traces
if cue == 1
    a = zeros(IntanBehaviour.nCueHit,1);
    shortTrialFlag = 0; % 0 - none, 1 - first, 2 = end , 3 = both
    
    for i=1:IntanBehaviour.nCueHit
        a(i) = max(find(cueIndex<=rewardIndex(i)));
        IntanBehaviour.cueHit(i,1) = cueIndex(a(i)); % Cue index 
        IntanBehaviour.cueHit(i,2) = rewardIndex(i); % pull index for hits
        if parameters.opto == 1
            if IntanBehaviour.optoTrace(IntanBehaviour.cueHit(i,1)) == 1   %~isempty(find(optoIndex == IntanBehaviour.cueHit(i,1))) 
                IntanBehaviour.cueHitTrace(i).opto = 1;
            else 
                IntanBehaviour.cueHitTrace(i).opto = 0;
            end 
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
    IntanBehaviour.reactionTime = (diff(IntanBehaviour.cueHit,1,2)/parameters.Fs)'; % in seconds
end

%% Getting cue Miss traces
if cue == 1
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
        if parameters.opto == 1
            if IntanBehaviour.optoTrace(IntanBehaviour.cueMiss(i,1)) == 1  %~isempty(find(optoIndex == IntanBehaviour.cueMiss(i,1))) 
                IntanBehaviour.cueMissTrace(i).opto = 1;
            else
                IntanBehaviour.cueMissTrace(i).opto = 0;
            end 
        end
        st = IntanBehaviour.cueMiss(i,1)-parameters.windowBeforeCue*parameters.Fs;
        sp = IntanBehaviour.cueMiss(i,1)+parameters.windowAfterCue*parameters.Fs;
        IntanBehaviour.cueMissTrace(i).trace = IntanBehaviour.leverTrace(st:sp)';
        IntanBehaviour.cueMissTrace(i).time = (0:1/parameters.Fs:(size(IntanBehaviour.cueMissTrace(i).trace,1)-1)*1/parameters.Fs)' - parameters.windowBeforeCue;
        IntanBehaviour.cueMissTrace(i).LFPIndex = ([st:1:sp])';
        IntanBehaviour.cueMissTrace(i).LFPtime = IntanBehaviour.time(st:sp)';
    end
end
%% Getting miss traces and timings

% figure();plot(1:1:10001,(5/1024)*Behaviour.leverTrace(Behaviour.miss(15,1)-5000:Behaviour.miss(15,1)+5000));hold on;
% plot(0.1:0.1:10000.1,IntanBehaviour.leverTrace(Behaviour.miss(15,3)-50000:Behaviour.miss(15,3)+50000));
correctionWindow = 1800; % in number of points in LFPFs
tol = 0.04;
nDiffSlope = 5;
disp('Finding miss trials in the Intan data ...');
for i=1:Behaviour.nMiss
    missIndexAr = Behaviour.miss(i,3);
    st1 = missIndexAr-correctionWindow;
    if st1 <= 0
        disp("Short miss trial removed");
        continue;
    end
    sp1 = missIndexAr+correctionWindow;
    if sp1 >= size(IntanBehaviour.leverTrace,2)
        disp("Short miss trial removed")
        continue;
    end
    trace1 = IntanBehaviour.leverTrace(st1:sp1);
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
    if ~isnan(misstrigs1)
        missIndex(i) =  missIndexAr - correctionWindow + min(misstrigs1);
    else
        missIndex(i) = NaN;
        disp('Correction failed. Trial removed');
    end
end

missIndex = removeNaNRows(missIndex');
idx =  missIndex<parameters.windowBeforePull*parameters.Fs;
missIndex(idx) = [];
IntanBehaviour.nMiss = size(missIndex,1);

for i=1:IntanBehaviour.nMiss
%     IntanBehaviour.hit(i) = [rewardIndex(i) lfpTime(rewardIndex(i)) rewardIndex(i) lfpTime(rewardIndex(i))];
    IntanBehaviour.missTrace(i).trace = IntanBehaviour.leverTrace(missIndex(i)-parameters.windowBeforePull*parameters.Fs:missIndex(i)+parameters.windowAfterPull*parameters.Fs)';
    IntanBehaviour.missTrace(i).time = (0:1/parameters.Fs:(size(IntanBehaviour.missTrace(i).trace,1)-1)*1/parameters.Fs)' - parameters.windowBeforePull;
    IntanBehaviour.missTrace(i).LFPIndex = ([missIndex(i)-parameters.windowBeforePull*parameters.Fs:1:missIndex(i)+parameters.windowAfterPull*parameters.Fs])';
    IntanBehaviour.missTrace(i).LFPtime = IntanBehaviour.time(missIndex(i)-parameters.windowBeforePull*parameters.Fs:missIndex(i)+parameters.windowAfterPull*parameters.Fs)';
end

%% Getting cue Hit traces with allignment at MI 
meanMotionTrace = mean(horzcat(IntanBehaviour.hitTrace.trace),2);
IntanBehaviour.meanRestingPositionCue = mean(meanMotionTrace(1:500));
IntanBehaviour.MIcutoffHit = 0.2*(IntanBehaviour.threshold-IntanBehaviour.meanRestingPositionCue) + IntanBehaviour.meanRestingPositionCue;

badTrials = [];

for i=1:size(IntanBehaviour.cueHit,1)
    f = IntanBehaviour.hitTrace(i).trace - IntanBehaviour.MIcutoffHit;
    fAbove = f.*(f >= 0);
    fCross = find(diff(fAbove>0)==1);
    rewardIndex = (parameters.windowBeforePull*parameters.Fs) + 1;
    fCross(fCross>rewardIndex) = [];
    try
    if isempty(fCross)
        badTrials = [badTrials,i];
        disp('Bad trial in Motion Allignment Detected');
%         IntanBehaviour.MIHitTrace(i).MIIndex = IntanBehaviour.hitTrace(i).LFPIndex(fCross(end));
%         IntanBehaviour.MIHitTrace(i).trace = IntanBehaviour.leverTrace(IntanBehaviour.MIHitTrace(i).MIIndex-parameters.windowBeforeMI*parameters.Fs:IntanBehaviour.MIHitTrace(i).MIIndex+parameters.windowAfterMI*parameters.Fs)';
%         IntanBehaviour.MIHitTrace(i).time = (0:1/parameters.Fs:(size(IntanBehaviour.MIHitTrace(i).trace,1)-1)*1/parameters.Fs)' - parameters.windowBeforeMI;
%         IntanBehaviour.MIHitTrace(i).LFPIndex = ([IntanBehaviour.MIHitTrace(i).MIIndex-parameters.windowBeforeMI*parameters.Fs:1:IntanBehaviour.MIHitTrace(i).MIIndex+parameters.windowAfterMI*parameters.Fs])';
%         IntanBehaviour.MIHitTrace(i).LFPtime = IntanBehaviour.time(IntanBehaviour.MIHitTrace(i).MIIndex-parameters.windowBeforeMI*parameters.Fs:IntanBehaviour.MIHitTrace(i).MIIndex+parameters.windowAfterMI*parameters.Fs)';
    else
        IntanBehaviour.MIHitTrace(i).MIIndex = IntanBehaviour.hitTrace(i).LFPIndex(fCross(end));
        IntanBehaviour.MIHitTrace(i).cueIndex = IntanBehaviour.cueHitTrace(i).LFPIndex(parameters.windowBeforeCue*parameters.Fs+1);
        IntanBehaviour.MIHitTrace(i).rewardIndex = IntanBehaviour.hitTrace(i).LFPIndex(parameters.windowBeforePull*parameters.Fs+1);
        IntanBehaviour.MIHitTrace(i).reactionTime = (1/parameters.Fs)*(IntanBehaviour.MIHitTrace(i).MIIndex-IntanBehaviour.MIHitTrace(i).cueIndex);
        IntanBehaviour.cueHitTrace(i).reactionTime = IntanBehaviour.MIHitTrace(i).reactionTime;
        IntanBehaviour.cueHitTrace(i).rewardIndex = IntanBehaviour.MIHitTrace(i).rewardIndex;
        IntanBehaviour.MIHitTrace(i).trace = IntanBehaviour.leverTrace(IntanBehaviour.MIHitTrace(i).MIIndex-parameters.windowBeforeMI*parameters.Fs:IntanBehaviour.MIHitTrace(i).MIIndex+parameters.windowAfterMI*parameters.Fs)';
        IntanBehaviour.MIHitTrace(i).time = (0:1/parameters.Fs:(size(IntanBehaviour.MIHitTrace(i).trace,1)-1)*1/parameters.Fs)' - parameters.windowBeforeMI;
        IntanBehaviour.MIHitTrace(i).LFPIndex = ([IntanBehaviour.MIHitTrace(i).MIIndex-parameters.windowBeforeMI*parameters.Fs:1:IntanBehaviour.MIHitTrace(i).MIIndex+parameters.windowAfterMI*parameters.Fs])';
        IntanBehaviour.MIHitTrace(i).LFPtime = IntanBehaviour.time(IntanBehaviour.MIHitTrace(i).MIIndex-parameters.windowBeforeMI*parameters.Fs:IntanBehaviour.MIHitTrace(i).MIIndex+parameters.windowAfterMI*parameters.Fs)';
    end
    catch ME
        badTrials = [badTrials,i];
        disp('Bad trial in Motion Allignment Detected');
        continue
    end
end

% Removing bad trials

IntanBehaviour.hitTrace(badTrials) = [];
IntanBehaviour.cueHitTrace(badTrials) = [];
IntanBehaviour.MIHitTrace(badTrials) = [];
IntanBehaviour.nHit = size(IntanBehaviour.hitTrace,2);
IntanBehaviour.nCueHit = size(IntanBehaviour.cueHitTrace,2);

%% Getting FA traces with allignment at MI

meanMotionTrace = mean(horzcat(IntanBehaviour.missTrace.trace),2);
IntanBehaviour.meanRestingPositionFA = mean(meanMotionTrace(1:500));
IntanBehaviour.MIcutoffFA = 0.2*(IntanBehaviour.threshold-IntanBehaviour.meanRestingPositionFA) + IntanBehaviour.meanRestingPositionFA;

for i=1:IntanBehaviour.nMiss
    f = IntanBehaviour.missTrace(i).trace - IntanBehaviour.MIcutoffFA;
    fAbove = f.*(f >= 0);
    fCross = find(diff(fAbove>0)==1);
    rewardIndex = (parameters.windowBeforePull*parameters.Fs) + 1;
    fCross(fCross>rewardIndex) = [];
    if isempty(fCross)
        fCross = rewardIndex;
    end
    IntanBehaviour. MIFATrace(i).MIIndex = IntanBehaviour.missTrace(i).LFPIndex(fCross(end));
    if IntanBehaviour. MIFATrace(i).MIIndex-parameters.windowBeforeMI*parameters.Fs<0
        continue;
    end
    IntanBehaviour. MIFATrace(i).trace = IntanBehaviour.leverTrace(IntanBehaviour. MIFATrace(i).MIIndex-parameters.windowBeforeMI*parameters.Fs:IntanBehaviour. MIFATrace(i).MIIndex+parameters.windowAfterMI*parameters.Fs)';
    IntanBehaviour. MIFATrace(i).time = (0:1/parameters.Fs:(size(IntanBehaviour. MIFATrace(i).trace,1)-1)*1/parameters.Fs)' - parameters.windowBeforeMI;
    IntanBehaviour. MIFATrace(i).LFPIndex = ([IntanBehaviour. MIFATrace(i).MIIndex-parameters.windowBeforeMI*parameters.Fs:1:IntanBehaviour. MIFATrace(i).MIIndex+parameters.windowAfterMI*parameters.Fs])';
    IntanBehaviour. MIFATrace(i).LFPtime = IntanBehaviour.time(IntanBehaviour. MIFATrace(i).MIIndex-parameters.windowBeforeMI*parameters.Fs:IntanBehaviour. MIFATrace(i).MIIndex+parameters.windowAfterMI*parameters.Fs)';
end
badTrials = arrayfun(@(x) isempty(x.trace), IntanBehaviour.MIFATrace);
IntanBehaviour.MIFATrace(badTrials) = [];
if plotOption == 1
    % Plotting Lever traces for Cue Hit 
    figure('Name','Average Lever Traces for Cue Hits and Misses');
    subplot(211);
    for i=1:size(IntanBehaviour.cueHitTrace,2)
        plot(IntanBehaviour.cueHitTrace(i).time,IntanBehaviour.cueHitTrace(i).trace,'Color',[0 0 0 0.1],'LineWidth',1.5);
        hold on;
    end
    plot(IntanBehaviour.cueHitTrace(1).time,mean(horzcat(IntanBehaviour.cueHitTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
    yline(IntanBehaviour.threshold,'--.b','Threshold','LabelHorizontalAlignment','left'); 
    xline(0,'--r','Cue','LabelVerticalAlignment','top');%ylim([0 0.1]);
    xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
    ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for Cue Hits');box off;
    
    % Plotting Lever traces for Cue Miss 
    subplot(212);
    for i=1:size(IntanBehaviour.cueMissTrace,2)
        plot(IntanBehaviour.cueMissTrace(i).time,IntanBehaviour.cueMissTrace(i).trace,'Color',[0 0 0 0.2],'LineWidth',1.5);
        hold on;
    end
    plot(IntanBehaviour.cueMissTrace(1).time,mean(horzcat(IntanBehaviour.cueMissTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
    yline(IntanBehaviour.threshold,'--.b','Threshold','LabelHorizontalAlignment','left'); 
    xline(0,'--r','Cue','LabelVerticalAlignment','top');
    ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for Cue Misses');box off;
    
    % Plotting Lever traces for Hits and False Alarms
    figure('Name','Average Lever Traces for Hits and False Alarms');
    subplot(2,1,1);
    for i=1:size(IntanBehaviour.hitTrace,2)
        plot(IntanBehaviour.hitTrace(i).time,IntanBehaviour.hitTrace(i).trace,'Color',[0 0 0 0.1],'LineWidth',1.5);
        hold on;
    end
    plot(IntanBehaviour.hitTrace(1).time,mean(horzcat(IntanBehaviour.hitTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
    yline(IntanBehaviour.threshold,'--.b','Threshold','LabelHorizontalAlignment','left'); 
    xline(0,'--r','Reward','LabelVerticalAlignment','top');%ylim([0 0.1]);
    ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for Hits');box off;
    
    subplot(2,1,2);
    for i=1:size(IntanBehaviour.missTrace,2)
        plot(IntanBehaviour.missTrace(i).time,IntanBehaviour.missTrace(i).trace,'Color',[0 0 0 0.1],'LineWidth',1.5);
        hold on;
    end
    plot(IntanBehaviour.missTrace(1).time,mean(horzcat(IntanBehaviour.missTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
    yline(IntanBehaviour.threshold,'--.b','Threshold','LabelHorizontalAlignment','left'); 
    xline(0,'--r','False Alarm','LabelVerticalAlignment','top');ylim([0 0.1]);
    ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for Hits');box off;
    
    % Plotting Lever traces for MI alligned cue Hit Trails 
    figure('Name','Average Lever Traces for Cue Hits and FAs with allignment at MI');
    subplot(2,1,1);
    for i=1:size(IntanBehaviour.MIHitTrace,2)
        plot(IntanBehaviour.MIHitTrace(i).time,IntanBehaviour.MIHitTrace(i).trace,'Color',[0 0 0 0.1],'LineWidth',1.5);
        hold on;
    end
    plot(IntanBehaviour.MIHitTrace(1).time,mean(horzcat(IntanBehaviour.MIHitTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
    yline(IntanBehaviour.MIcutoffHit,'--.b','MI Threshold','LabelHorizontalAlignment','left'); 
    xline(0,'--r','MI','LabelVerticalAlignment','top');ylim([0 0.1]);
    % xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
    ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for Cue Hits');box off;
    
    subplot(2,1,2);
    for i=1:size(IntanBehaviour.MIFATrace,2)
        plot(IntanBehaviour.MIFATrace(i).time,IntanBehaviour.MIFATrace(i).trace,'Color',[0 0 0 0.1],'LineWidth',1.5);
        hold on;
    end
    plot(IntanBehaviour.MIFATrace(1).time,mean(horzcat(IntanBehaviour.MIFATrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
    yline(IntanBehaviour.MIcutoffFA,'--.b','MI Threshold','LabelHorizontalAlignment','left'); 
    xline(0,'--r','MI','LabelVerticalAlignment','top');ylim([0 0.1]);
    % xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
    ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for False Alarms');box off;
end
