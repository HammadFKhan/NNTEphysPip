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
channel_num = 1:size(amplifier_data,1);
% Notch filtering at 60 Hz
Fs = 20000;
d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',Fs);
% Design butterworth filters for single unit
Fc = [250 1000];
Wn = Fc./(Fs/2);
[b1,a1] = butter(6,Wn,'bandpass');

% Design butterworth filters for LFP
Fc = [250];
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
%%
intanData = preprocess_filtering(amplifier_data);
intanData = spikeSorting(intanData);
%%
for i = 1:channel_num(end)
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
fn = fieldnames(dataLow);
% findData = zeros(30,.3*Fs+1);
count = 1;
for i = 1:size(rippleBatch,2)
    if size(fieldnames(rippleBatch(count).rippleData),1) < 3
        rippleBatch(count) = []; 
    else
        count = count+1;
    end
end

for ii = 1:size(rippleBatch,2)
    LFPData = dataLow.(string(fn(ii)));
    rawData = intanData.rawData(ii,:);
    findCh = (rippleBatch(ii).rippleData.ripples(2,2))*Fs;
    try
        findDataLFP(ii,:) = LFPData(1,-.1*Fs+findCh:.1*Fs+findCh);
        findDataRaw(ii,:) = rawData(1,-.1*Fs+findCh:.1*Fs+findCh);
    catch ME
        disp('Not all channels were analyzed')
        break
    end
end
figure('Name','SWR Onset Frequency')
for j = 1:size(findDataRaw,1)
    [cfs,f] = cwt(findDataRaw(j,:),Fs,'FrequencyLimits',[20 300]);
    % plot(abs(cfs(1,:)));
    subplot(6,6,j),imagesc(-100:100,f,abs(cfs))
%     xlabel('Time (s)')
%     ylabel('Frequency (Hz)')
    axis xy
    colormap(jet)
    % ylim([0 250])
%     title('CWT of Ripple Data')
end
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
    print(gcf,'-painters','-depsc', 'Figures/SWRonset.eps', '-r250');

figure('Name','CSD');
channelmap = ones(size(findDataLFP,1),1);
try
    Vq = interp2(findDataLFP,5);
catch ME
    disp('Sample set is too small for interpolation, plotting raw...');
    Vq = findDataLFP;
end
imagesc(-100:100,channelmap,Vq);colormap(jet); colorbar;box off;set(gca,'YTick',[]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
    print(gcf,'-painters','-depsc', 'Figures/CSD.eps', '-r250');
%% Plots
disp('Plotting...')
figure('Name', 'Unfiltered Data'),stack_plot(amplifier_data);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
    print(gcf,'-painters','-depsc', 'Figures/Raw_data.eps', '-r250');
figure('Name','Singe Unit Waveforms'),SingleUnits(allSpikes);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
    print(gcf,'-painters','-depsc', 'Figures/Single_unit_waveforms.eps', '-r250');
figure('Name','Multi-Unit Activity'),MUA(dataHigh);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
    print(gcf,'-painters','-depsc', 'Figures/MultiUnit.eps', '-r250');
figure('Name','LFP'),LFP(dataLow); 
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
    print(gcf,'-painters','-depsc', 'Figures/LFP.eps', '-r250');
%%
try
encoder_data = convert_encoder(board_adc_data(2,:),timestamps);
catch ME
    disp('No ADC data')
    return
end
figure('Name','Pulse Data');plot(encoder_data.rotate_pulse);
figure('Name','Angular Distance');bar(encoder_data.ang_distance);
figure('Name','Angular Velocity');bar(encoder_data.ang_velocity,'FaceColor',[.16 .835 .384],'EdgeColor','none');
figure('Name','Avg. Angular Velocity');avgV = movmean(encoder_data.ang_velocity,2);bar(avgV,'FaceColor',[.16 .835 .384],'EdgeColor','none');
