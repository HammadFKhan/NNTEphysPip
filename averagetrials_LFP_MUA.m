clear
%get good channels
channels = 1:6;
Fs = 20000;

%get names of all intan files
fileList = dir('*.rhd');
numTrials = length(fileList);
% numTrials = 20;
path = pwd;

%define window around pole movement
windowBefore = .025*Fs;
windowAfter = 0.250*Fs;

numStd = 5;

%%
%design butterworth filters for single unit
Fc = [300 3000];
Wn = Fc./(Fs/2);
[b1,a1] = butter(6,Wn,'bandpass');

%design butterworth filters for LFP
Fc = [250];
Wn = Fc./(Fs/2);
[b2,a2] = butter(6,Wn,'low');



%%
count = 0;
spikeCount = zeros(length(channels),1);
for tt = 1:numTrials
    %load the first intan file
    read_Intan_RHD2000_fileName([path '\'],fileList(tt).name)
    
    %find beginning of pole movement
    poleMoveIndices = find(board_dig_in_data(2,:));
    dd = find(diff(board_dig_in_data(2,:)));
    touchStartIndex = dd(1);
    poleLength = length(dd);
    
    trialIndices = touchStartIndex-windowBefore:touchStartIndex+windowAfter;
    trialTime = (1:length(trialIndices))./Fs.*1000;
    poleTime   = (windowBefore:(windowBefore+length(poleMoveIndices)))./Fs.*1000;
    
    for cc = 1:length(channels)
        count = count+1;
        
%         %notch filer
%         notchData = filter(b3,a3,amplifier_data(channels(cc),:));
%         notchData = filter(b4,a4,notchData,:);

        %to get LFP band, first highpass the data, then lowpass
        %highpass data to get rid of really low freq stuff (<10Hz)
        [highpassData,d1] = highpass(amplifier_data(channels(cc),:),5,Fs,'ImpulseResponse','iir','Steepness',0.5);
        lowpassData  = filtfilt(b2,a2,double(highpassData));
        
        %save the LFP trace
        dataLow.(['Channel' num2str(cc)])(tt,:) = lowpassData(trialIndices);
        
        %bandpass data for spike band
        bandpassData = filtfilt(b1,a1,double(amplifier_data(channels(cc),:)));
        
        %save bandpassed data for plotting
        dataHigh.(['Channel' num2str(cc)])(tt,:) = bandpassData(trialIndices);
        
        %set a threshold
        thresh = numStd.*std(bandpassData);
        %only look at the bandpassed data during the trial start/stop
        bandpassData = bandpassData(trialIndices);
        
        %find spikes
        [pks,locs] = findpeaks(-bandpassData,Fs,'MinPeakHeight',thresh,'MinPeakDistance',.001);
        
        %for each detected spike, extract a window 2ms before and after the peak
        window = .00175*Fs;%5ms
        if(~isempty(pks))
            allLocs.(['Channel' num2str(channels(cc))]){tt} = locs;
            for pp = 1:length(locs)
                thisLoc = locs(pp)*Fs;
                if(thisLoc-window>1 && thisLoc+window<length(bandpassData))%%peak can't occur at the very end or very beginning of the data set
                   spikeCount(channels(cc)) = spikeCount(channels(cc))+1;
                    %extract a 4ms window around the spike peak
                    allSpikes.(['Channel' num2str(channels(cc))])(spikeCount(channels(cc)),:) = bandpassData(thisLoc-window:thisLoc+window);        
                end
            end
        else
            allLocs.(['Channel' num2str(channels(cc))]){tt} = [];
        end
    end
end


%% plot LFP
spacing = 1000;
ycenter =[spacing:spacing:spacing*length(channels)]; 
colorgrad = zeros(length(channels),3);
colorgrad(:,1) = linspace(0.6350,0     ,length(channels));
colorgrad(:,2) = linspace(0.0780,0.4470,length(channels));
colorgrad(:,3) = linspace(0.1840,0.7410,length(channels));

figure;hold on

for cc = 1:length(channels)
    %plot LFP
    channelData = dataLow.(['Channel' num2str(cc)]);
    plot(trialTime,channelData'+ycenter(cc),'color',colorgrad(cc,:))
    plot(trialTime,mean(channelData)+ycenter(cc),'color','k','LineWidth',2)
    
%     %slightly under it, plot MUA
%     channelData = dataHigh.(['Channel' num2str(cc)]);
%     plot(trialTime,channelData'+ycenter(cc)-200,'color',colorgrad(cc,:))
%     plot(trialTime,mean(channelData)+ycenter(cc)-200,'color','k','LineWidth',2)    
end

rectangle('Position',[windowBefore./Fs.*1000,spacing,length(poleMoveIndices)/Fs*1000,spacing*length(channels)],...
          'FaceColor',[0.25 0.25 0.25 0.2],...
          'EdgeColor','none')
alpha(0.2)

box off 
set(gca,'YTick','','YTickLabel','')
title('LFP Freq Band')
xlabel('Time (ms)')
ylabel('Channels')


%% plot MUA
spacing = 100;
ycenter =[spacing:spacing:spacing*length(channels)]; 
colorgrad = zeros(length(channels),3);
colorgrad(:,1) = linspace(0.6350,0     ,length(channels));
colorgrad(:,2) = linspace(0.0780,0.4470,length(channels));
colorgrad(:,3) = linspace(0.1840,0.7410,length(channels));


figure;hold on

for cc = 1:length(channels)
    channelData = dataHigh.(['Channel' num2str(cc)]);

    plot(trialTime,channelData'+ycenter(cc),'color',colorgrad(cc,:))
    plot(trialTime,mean(channelData)+ycenter(cc),'color','k','LineWidth',2)
%         if(tt==numTrials)
%             plot(poleTime,ycenter(cc).*ones(1,length(poleTime)),'color',[0.25, 0.25, 0.25],'LineWidth',5)
%         end
end

rectangle('Position',[windowBefore./Fs.*1000,spacing,length(poleMoveIndices)/Fs*1000,spacing*length(channels)],...
          'FaceColor',[0.25 0.25 0.25 0.2],...
          'EdgeColor','none')
alpha(0.2)

box off 
set(gca,'YTick','','YTickLabel','')
title('LFP Freq Band')
xlabel('Time (ms)')
ylabel('Channels')

%% plot single unit waveforms
figure
for cc = 1:length(channels)
    subplot(1,length(channels),cc);hold on
    theseSpikes = allSpikes.(['Channel' num2str(channels(cc))]);
    
    t = (1:length(theseSpikes(1,:)))./Fs.*1000;
    plot(t,theseSpikes,'b')
    hold on
    plot(t,mean(theseSpikes),'k','LineWidth',2)
end
   
