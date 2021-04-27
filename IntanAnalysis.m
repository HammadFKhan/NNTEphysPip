%% Intan Data
clear; clc; close all;
read_Intan_RHD2000_file
% amplifier_data = amplifier_data(:,1:100000);
% t_amplifier = t_amplifier(1:100000);
%%
channel_num = 15:31;
figure()
plot(t_amplifier,amplifier_data(channel_num,:));
%%
lowpassfilter = [lowpass(amplifier_data(channel_num,:),10,frequency_parameters.amplifier_sample_rate/5,'ImpulseResponse','iir','Steepness',0.75)];
figure()
plot(t_amplifier,lowpassfilter);
%%
% highpassfilter = [highpass(amplifier_data(channel_num,:),500,frequency_parameters.amplifier_sample_rate,'ImpulseResponse','iir','Steepness',0.50)];
% figure()
% plot(t_amplifier,highpassfilter);
% %%
% bandpassfilter = [bandpass(amplifier_data(channel_num,:),[300,3000],frequency_parameters.amplifier_sample_rate,'ImpulseResponse','iir','Steepness',0.60)];
% figure()
% plot(t_amplifier,bandpassfilter);
%% Butterworth
for i = 1:10
    rawData = amplifier_data(channel_num(i),:);
    Fs = 20000;
    Fc = [300 3000];
    Wn = Fc./(Fs/2);
    [b1,a1] = butter(4,Wn,'bandpass');
    bandpassData(i,:) = filtfilt(b1,a1,double(rawData));
    figure(i)
    plot(t_amplifier,bandpassData(i,:));
    % Spike Detection
    dt = t_amplifier(2)-t_amplifier(1);
    ds(i,:) = diff(bandpassData(i,:));
    ds_sub(i) = mean(ds(i,floor(0.5/dt)+1:floor(1/dt)));
    ds_std(i) = std(ds(i,floor(0.5/dt)+1:floor(1/dt)));
    threshold(i) = 3*ds_std(i);
    spike_detection_ds{i} = [];
    for n = 2:length(ds)
        if abs(ds(i,n))>threshold(i) && abs(ds(i,n-1))<=threshold(i)
            spike_detection_ds{i} = [spike_detection_ds{i}, n*dt]; % spike time in s
%         else 
%             spike_detection_ds{i} = [spike_detection_ds{i}, NaN];
        end
    end
end

%% Plot Spike Plot
% tic
% figure('Name','Spike Plot','NumberTitle','off')
% for x = 1:length(spike_detection_ds)
%     dspike_vector(x,1) = length(spike_detection_ds{x});
% end
% dspike_vector_max = max(dspike_vector);
%     for j = 1:length(spike_detection_ds)
%         if length(spike_detection_ds{j}) ~= dspike_vector_max
%             dif = abs(length(spike_detection_ds{j}) - dspike_vector_max);
%             spike_detection_ds{j}(end+dif) = 0;
%         end
%         for n = 1:dspike_vector_max
%             if spike_detection_ds{j}(n) == 0
%                 spike_detection_ds{j}(n) = NaN;
%             end
%         end
%         plot(spike_detection_ds{j},j,'.','MarkerSize',10,'Color','black'); hold on;
%         dspike_plot(:,j) = cell2mat(spike_detection_ds(j)'); 
%         
%     end
%     axis([0 t_amplifier(end) 0 length(spike_detection_ds)+1]);
%     
% % figure('Name','Raster Plot','NumberTitle','off')
% %  for v = 1:length(spike_detection_ds)
% %      plot_raster(spike_detection_ds{v},v); hold on;
% %  end
% %  xlim([0 T(end)]);
% %  ylim([0 length(spike_detection_ds)+2]);
% %  ax = gca;
%  ax.Box = 'off';
%  rastertime = toc
% 
% 
%% find all spikes and align to the peak
%find negative peaks in filtered data with an amplitude of at least 4 std and 1ms  in
%separation
filteredData = lowpassfilter;
thresh = 3*std(filteredData);
[pks,locs] = findpeaks(-filteredData,Fs,'MinPeakHeight',thresh,'MinPeakDistance',.001);
%for each detected spike, extract a window 2ms before and after the peak
spikeCount = 0;
window = .00175*Fs;%5ms
if(~isempty(pks))
    for pp = 1:length(locs)
        thisLoc = locs(pp)*Fs;
        if(thisLoc-window>1 && thisLoc+window<length(filteredData))%%peak can't occur at the very end or very beginning of the data set
            spikeCount = spikeCount+1;
            %extract a 4ms window around the spike peak
            allSpikes(spikeCount,:) = filteredData(thisLoc-window:thisLoc+window);
        end
    end
end
figure()
plot(all_Spikes);