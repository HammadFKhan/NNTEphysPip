%% Intan Data
clear; clc; 
close all;
read_Intan_RHD2000_file
% % amplifier_data = amplifier_data(:,1:100000);
% t_amplifier = t_amplifier(1:100000);
%%
channel_num = 1:31;
% Design butterworth filters for single unit
Fs = 20000;
Fc = [250 3000];
Wn = Fc./(Fs/2);
[b1,a1] = butter(6,Wn,'bandpass');

% Design butterworth filters for LFP
Fc = [25];
Wn = Fc./(Fs/2);
[b2,a2] = butter(6,Wn,'low');
% Design butterworth filter for Ripples
Fc = [180 220];
Wn = Fc./(Fs/2);
[b3,a3]=butter(4,Wn,'bandpass');
%% General Filtering
% Fc = [3000];
% Wn = Fc./(Fs/2);
% [b,a] = cheby2(4,20,Wn,'low');
% amplifier_data = filtfilt(b,a,amplifier_data);


%%
findSpikes = true;
spikeWindow = .00025*Fs;%35ms
threshStd = 5;
spikeCount = zeros(length(channel_num),1);
%%
H = waitbar(0,'Compiling...');
for i = 1:channel_num(end)
    waitbar(i/channel_num(end),H)
    rawData = amplifier_data(channel_num(i),:);
    lowpassData  = filtfilt(b2,a2,rawData);
    bandpassData = filtfilt(b1,a1,rawData);
    
    data.rawData(i,:) = rawData;
    data.lowpassData(i,:) = lowpassData;
    data.bandpassData(i,:) = bandpassData;
    
    dataLow.(['Channel' num2str(i)]) = lowpassData;
    dataHigh.(['Channel' num2str(i)]) = bandpassData;
    
    %find spikes
    thresh = threshStd*std(bandpassData);
    [pks,locs] = findpeaks(-bandpassData,Fs,'MinPeakHeight',thresh,'MinPeakDistance',.001);
    %save the spike times for this channel
    window = .00175*Fs;%35ms
    if(~isempty(pks))
        allLocs.(['Channel' num2str(channel_num(i))]) = locs;
        for pp = 1:length(locs)
            thisLoc = locs(pp)*Fs;
            if(thisLoc-window>1 && thisLoc+window<length(bandpassData))%%peak can't occur at the very end or very beginning of the data set
                spikeCount(channel_num(i)) = spikeCount(channel_num(i))+1;
                %extract a 4ms window around the spike peak
                allSpikes.(['Channel' num2str(channel_num(i))])(spikeCount(channel_num(i)),:) = bandpassData(thisLoc-window:thisLoc+window);
            end
        end
    else
        allLocs.(['Channel' num2str(channel_num(i))]) = [];
    end
    
    % Ripples
    ripple_signal(:,i) = FiltFiltM(b3,a3,rawData);
    pow(:,i) = fastrms(ripple_signal(:,i),15);
    mRipple(i) = mean(pow(:,i));
    meRipple(i) = median(pow(:,i));
    mmRippleRatio(i) = mRipple(i)./meRipple(i);
    % Ripple detection
    ripples =  rippleDetection(ripple_signal(:,i),Fs);
    allRipples{i} = ripples;
end
close(H)
    
mmRippleRatio(mRipple<1) = 0;
mmRippleRatio(meRipple<1) = 0;

[minVal,loc] = max(mmRippleRatio);
chan = channel_num(loc);
disp([ 'Best channel detected: ' num2str(chan)]);
%% Plots
disp('Plotting...')
figure('Name', 'Unfiltered Data'),stack_plot(amplifier_data)
figure('Name','Singe Unit Waveforms'),SingleUnits(allSpikes,Fs,125);
figure('Name','Multi-Unit Activity'),MUA(dataHigh);
figure('Name','LFP'),LFP(dataLow);