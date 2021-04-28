%% Intan Data
clear; clc; close all;
read_Intan_RHD2000_file
% amplifier_data = amplifier_data(:,1:100000);
% t_amplifier = t_amplifier(1:100000);
%%
channel_num = 1:31;
figure()
plot(t_amplifier,amplifier_data(channel_num,:));
%%
%design butterworth filters for single unit
Fs = 20000;
Fc = [250 3000];
Wn = Fc./(Fs/2);
[b1,a1] = butter(6,Wn,'bandpass');

%design butterworth filters for LFP
Fc = [50];
Wn = Fc./(Fs/2);
[b2,a2] = butter(6,Wn,'low');
%%
findSpikes = true;
spikeWindow = .00175*Fs;%5ms
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
    window = .00175*Fs;%5ms
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
end
close(H)

%% plot single unit waveforms

fn = fieldnames(allSpikes);
for ii = 1:11
%     subplot(1,20,ii);hold on
    theseSpikes = allSpikes.(string(fn(ii)));
    t = (1:length(theseSpikes(1,:)))./Fs.*1000;
    plot(t,theseSpikes,'b')
    hold on
    plot(t,mean(theseSpikes),'k','LineWidth',2)
end

%% plot LFP
time = (1:length(LFPData))./Fs.*1000;
fn = fieldnames(dataLow);
spacing = 1000;
ycenter =[spacing:spacing:spacing*length(channel_num)]; 
colorgrad = zeros(length(channel_num),3);
colorgrad(:,1) = linspace(0.6350,0     ,length(channel_num));
colorgrad(:,2) = linspace(0.0780,0.4470,length(channel_num));
colorgrad(:,3) = linspace(0.1840,0.7410,length(channel_num));

figure;hold on

for j = 14:length(channel_num)
    %plot LFP
    LFPData = dataLow.(string(fn(j)));
    plot(LFPData'+ycenter(j),'color',colorgrad(j,:))
    plot(mean(LFPData)+ycenter(j),'color','k','LineWidth',2)
    
%     %slightly under it, plot MUA
%     channelData = dataHigh.(['Channel' num2str(cc)]);
%     plot(trialTime,channelData'+ycenter(cc)-200,'color',colorgrad(cc,:))
%     plot(trialTime,mean(channelData)+ycenter(cc)-200,'color','k','LineWidth',2)    
end
% 
% rectangle('Position',[windowBefore./Fs.*1000,spacing,length(poleMoveIndices)/Fs*1000,spacing*length(channels)],...
%           'FaceColor',[0.25 0.25 0.25 0.2],...
%           'EdgeColor','none')
% alpha(0.2)

box off 
set(gca,'YTick','','YTickLabel','')
title('LFP Freq Band')
xlabel('Time (ms)')
ylabel('Channels')
xlim([3E5 5E5]);
%% plot MUA
time = (1:length(channelData))./Fs.*1000;
fn = fieldnames(dataHigh);
spacing = 100;
ycenter =[spacing:spacing:spacing*length(channel_num)]; 
colorgrad = zeros(length(channel_num),3);
colorgrad(:,1) = linspace(0.6350,0     ,length(channel_num));
colorgrad(:,2) = linspace(0.0780,0.4470,length(channel_num));
colorgrad(:,3) = linspace(0.1840,0.7410,length(channel_num));


figure;hold on

for jj = 14:length(channel_num)
    channelData = dataHigh.(string(fn(jj)));

    plot(channelData'+ycenter(jj),'color',colorgrad(jj,:))
    plot(mean(channelData)+ycenter(jj),'color','k','LineWidth',2)
%         if(tt==numTrials)
%             plot(poleTime,ycenter(cc).*ones(1,length(poleTime)),'color',[0.25, 0.25, 0.25],'LineWidth',5)
%         end
end

% rectangle('Position',[windowBefore./Fs.*1000,spacing,length(poleMoveIndices)/Fs*1000,spacing*length(channels)],...
%           'FaceColor',[0.25 0.25 0.25 0.2],...
%           'EdgeColor','none')
% alpha(0.2)

box off 
set(gca,'YTick','','YTickLabel','')
title('MUA Activity')
xlabel('Time (ms)')
ylabel('Channels')
xlim([3E5 5E5]);