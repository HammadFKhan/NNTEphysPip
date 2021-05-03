%% Intan Data
clear; clc; 
close all;
read_Intan_RHD2000_file
% % amplifier_data = amplifier_data(:,1:100000);
% t_amplifier = t_amplifier(1:100000);
amplifier_data_sorted = channelSortEdge(amplifier_data);
amplifier_data = amplifier_data_sorted;
%%
% set(0,'DefaultFigureWindowStyle','docked')
channel_num = 1:31;
% Design butterworth filters for single unit
Fs = 20000;
Fc = [250 3000];
Wn = Fc./(Fs/2);
[b1,a1] = butter(6,Wn,'bandpass');

% Design butterworth filters for LFP
Fc = [1];
Wn = Fc./(Fs/2);
[b2,a2] = butter(1,Wn,'high');
Fc = [50];
Wn = Fc./(Fs/2);
[b4,a4] = butter(2,Wn,'low');
% Design butterworth filter for Ripples
Fc = [140 200];
Wn = Fc./(Fs/2);
[b3,a3]=butter(4,Wn,'bandpass');
%% General Filtering
% Fc = [3000];
% Wn = Fc./(Fs/2);
% [b,a] = cheby2(4,20,Wn,'low');
% amplifier_data = filtfilt(b,a,amplifier_data);


%%
findSpikes = true;
spikeWindow = .00025*Fs;%5ms
threshStd = 5;
spikeCount = zeros(length(channel_num),1);
ripples = {};
%%
H = waitbar(0,'Compiling...');
for i = 1:channel_num(end)
    warning('off','all')
    waitbar(i/channel_num(end),H)
    rawData = amplifier_data(channel_num(i),:);
    highpass_LFP = filtfilt(b2,a2,rawData);
    lowpassData  = filtfilt(b4,a4,highpass_LFP);
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
    timestamps(:,1) = t_amplifier;
    ripple_signal(:,i) = FiltFiltM(b3,a3,rawData);
    pow(:,i) = fastrms(ripple_signal(:,i),15);
    mRipple(i) = mean(pow(:,i));
    meRipple(i) = median(pow(:,i));
    mmRippleRatio(i) = mRipple(i)./meRipple(i);
    % Ripple detection
    %     [ripples] =  rippleDetection(ripple_signal(:,i),timestamps,Fs)
    %     allRipples{i} = ripples;
    
    [rippleData,data,stats,maps] =  rippleDetection(ripple_signal(:,i),timestamps,Fs);
    rippleBatch(i).rippleData = rippleData;
    rippleBatch(i).data = data;
    rippleBatch(i).stats = stats;
    rippleBatch(i).maps = maps;
end
close(H)
mmRippleRatio(mRipple<1) = 0;
mmRippleRatio(meRipple<1) = 0;

[minVal,loc] = max(mmRippleRatio);
chan = channel_num(loc);
disp([ 'Best channel detected: ' num2str(chan)]);
% ripples = {};
ripplePlot(rippleData,data,stats,maps);
%%
figure()
[cfs,f] = cwt(amplifier_data(4,:)',Fs,'FrequencyLimit',[1 300]);
% plot(abs(cfs(1,:)));
imagesc(-100:100,f,abs(cfs))
xlabel('Time (s)')
ylabel('Frequency (Hz)')
axis xy
colormap(jet)
% ylim([0 250])
title('CWT of Ripple Data')
%%
fn = fieldnames(dataLow);
findData = zeros(30,.3*Fs+1);
count = 1;
for i = 1:size(rippleBatch,2)
    if size(fieldnames(rippleBatch(count).rippleData),1) < 1
        rippleBatch(count) = []; 
        count = count-1;
    else
        count = count+1;
    end
end

for ii = 1:size(rippleBatch,2)
    LFPData = dataLow.(string(fn(ii)));
    findCh = (rippleBatch(ii).rippleData.ripples(3,2))*Fs;
    findData(ii,:) = LFPData(1,-.150*Fs+findCh:.15*Fs+findCh);
end
Vq = interp2(findData,5);
figure;
imagesc(Vq);colormap(jet);
%% Plots
disp('Plotting...')
figure('Name', 'Unfiltered Data'),stack_plot(amplifier_data)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
    print(gcf,'-painters','-depsc', 'Figures/Raw_data.eps', '-r250');
figure('Name','Singe Unit Waveforms'),SingleUnits(allSpikes,Fs,200);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
    print(gcf,'-painters','-depsc', 'Figures/Single_unit_waveforms.eps', '-r250');
figure('Name','Multi-Unit Activity'),MUA(dataHigh);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
    print(gcf,'-painters','-depsc', 'Figures/MultiUnit.eps', '-r250');
figure('Name','LFP'),LFP(dataLow);xlim([5E5 7E5]); 
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
    print(gcf,'-painters','-depsc', 'Figures/LFP.eps', '-r250');