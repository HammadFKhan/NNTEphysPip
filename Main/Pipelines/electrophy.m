close all
clear all
idx = 2;

%  path = 'Z:\data\Hammad\01052021\';
%  filename = sprintf('210105_001.11.wcp', idx);
% path = 'Z:\data\Hammad\12142020\';
% filename = sprintf('201229_001.sub and supra.2.wcp', idx);
% path = 'Z:\data\Hammad\12292020\';
% filename = sprintf('201229_001.1.wcp', idx);
% path = 'Z:\data\Hammad\01162021\Cell2\';
% filename = sprintf('210116_001.LFPstim.8.wcp', idx);
% path = 'Z:\data\Hammad\01222021\spontaneous_1\';
% filename = sprintf('210122_001.2.wcp', idx);
% path = 'Z:\data\Hammad\01252021\spontaneous_sp\';
% filename = sprintf('210125_001.7.wcp', idx);
%% Sine Input
% path = 'Z:\data\Hammad\12122020\CA1_2';
% filename = sprintf('201212_001.Sine.2.wcp', idx);
% filename = sprintf('201214_001.Sine.2.wcp', idx);
% path = 'Z:\data\Hammad\12162020\';
% filename = sprintf('201216_001.Sine.1.wcp', idx);
% path = 'Z:\data\Hammad\12292020\';
tic
% out=import_wcp(fullfile(path, filename),'debug');
out=import_wcp;
n_channel = out.channel_no;
n_recording = length(out.rec_index);
% n_recording = 1;
%% 
dt = out.T(2)-out.T(1); % sample time in s
% spike_detection = {}; % spike time sorted as the rising phase when cross threshold
thr = 0; % voltage threshold for spike detection, in mV
H1 = waitbar(0,'Loading Ephys Data');
cutoff = 100;
std_val = 5;
for i = 8
    waitbar(i/n_recording)
    Vm(:,i) = out.S{3}(:,i);
    Im(:,i) = out.S{4}(:,i);
%     spike_detection{i} = [];
%     for n = floor(2/dt)+1:floor(2.5/dt)
%         if Vm(n, i)>thr && Vm(n-1, i)<=thr
%             spike_detection{i} = [spike_detection{i}, n*dt]; % spike time in s
%         end
%     end
    Ihold = mean(Im(1:floor(0.5/dt),i)); % holding current, in pA
    Istep_sub(i) = floor((mean(Im(floor(0.5/dt)+1:floor(1/dt),i)) - Ihold)/10)*10; % injected step current, in pA
    Istep_supra(i) = floor((mean(Im(floor(2/dt)+1:floor(2.5/dt),i)) - Ihold)/50)*50; % injected step current, in pA
    V_sub(i) = mean(Vm(floor(0.5/dt)+1:floor(1/dt),i));
    Vrest(i) = mean(Vm(1:floor(0.5/dt) ,i));
    if Istep_sub(i)== 0
        Rin(i) = 0;
    else
        Rin(i) = (V_sub(i)-Vrest(i))/Istep_sub(i)*1000; % input resistance, MOhm
    end
%     FR(i) = length(spike_detection{i})/0.5; % firing rate, Hz
    
    % lowpass filtering
%     [lowpass_Vm_plot{i},d] = [lowpass(Vm(:,i), 500,(1/(out.t_interval)),'ImpulseResponse','iir','Steepness',0.95)]; 
    lowpass_Vm(:,i) = [lowpass(Vm(:,i), cutoff,(1/(out.t_interval)),'ImpulseResponse','iir','Steepness',0.75)];
%   subplot(3,1,3),periodogram(lowpass_Vm(:,N));
    % Fourier Recovery
    L = length(Vm(:,i));
    Fs = 1/out.t_interval;
    fourier = fft(Vm(:,i));
    P2 = abs(fourier/L);
    P1 = (1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    % Power Spectrum
%     psdest = psd(spectrum.periodogram,Vm(:,N),'Fs',Fs,'NFFT',L);
%     plot(psdest.Frequencies, psdest.Data); grid on;
    
    % Derivative
    ds(:,i) = diff(Vm(:,i));
    secondds(:,i) = diff(Vm(:,i),2);
    ds_lowpass(:,i) = diff(lowpass_Vm(:,i));
    secondds_lowpass(:,i) = diff(lowpass_Vm(:,i),2);
%     Spikes = spike_detection(Vm,dt,std_val);
%     Spike detection using derivative
    ds_sub(i) = mean(ds_lowpass(floor(0.5/dt)+1:floor(1/dt),i));
    ds_std(i) = std(ds_lowpass(floor(0.5/dt)+1:floor(1/dt),i));
    threshold(i) = std_val*ds_std(i);
    spike_detection_ds{i} = [];
    for n = 2:length(ds_lowpass)
        if abs(ds_lowpass(n,i))>threshold(i) && abs(ds_lowpass(n-1,i))<=threshold(i)
            spike_detection_ds{i} = [spike_detection_ds{i}, n*dt]; % spike time in s
%         else 
%             spike_detection_ds{i} = [spike_detection_ds{i}, NaN];
        end
    end
   
%     spike_detection_ds = Spikes;
    dFR(i) = length(spike_detection_ds{i})/0.5;
end
Ihold = sum(mean(Im(1:floor(0.5/dt),:))/n_recording); % holding current, in pA
Vrest = sum(mean(Vm(1:floor(0.5/dt),:))/n_recording);
Rin_pooled = mean(Rin(find(Istep_sub~=0))); % input resistance, MOhm
T= out.T';
delete(H1)

% Plot all
figure()
for N = 8
%     figure('Name',strcat('Subplots',int2str(N)),'NumberTitle','off');dplot(T,Vm(:,N),lowpass_Vm(:,N),ds(:,N),ds_lowpass(:,N),...
%     secondds(:,N),secondds_lowpass(:,N));
% 
% %    Im_Plt = subplot(2,1,2), plot(T, Im(:,N), 'Linewidth',1); axis off;
% %    saveas(gcf,strcat('01052021Analysis\',filename,num2str(N),'.tif'));
% %    print(gcf,'uEPSPt','-dtiff','-r600');
% end

   % Wavelet Transform
    cwt(Vm(:,N),Fs);
%     [WT,F] = cwt(Vm(:,ii),Fs);
%     xrec = icwt(WT,F,[40,150],'SignalMean',mean(Vm(:,ii)));
%     plot(T,xrec);
%     ylim([floor(min(Vm(:,ii))) ceil(max(Vm(:,ii)))]);
%     xlim([0 T(end)]);

end
toc
figure()
for nn = 8
subplot(2,1,1),plot(T,lowpass_Vm(:,nn)); hold on;
ax = gca;
ax.Box = 'off';
subplot(2,1,2),plot(T(2:end),ds_lowpass(:,nn));hold on;
ax = gca;
ax.Box = 'off';
xlim([0 T(end)]);
end
% figure()
% for nnn = 3:4
% plot(T(2:end),ds_lowpass(:,nnn));hold on;
% ax = gca;
% ax.Box = 'off';
% xlim([0 T(end)]);
% end
% figure()
% for ii = 1:n_recording
%     highpassfilter(:,ii) = [highpass(lowpass_Vm(:,ii),500,1/out.t_interval,'ImpulseResponse','iir','Steepness',0.50)]; 
%     plot(T,highpassfilter(:,ii));hold on;
% end

for ii = 8
    rawData = Vm(:,ii);
    Fs = (1/(out.t_interval));
    Fc = [300 3000];
    Wn = Fc./(Fs/2);
    [b1,a1] = butter(4,Wn,'bandpass');
    bandpassData(ii,:) = filtfilt(b1,a1,double(rawData));
    figure()
    subplot(2,1,1),plot(T,bandpassData(ii,:));
end

tic
%Spike Plot
figure('Name','Spike Plot','NumberTitle','off')
for x = 1:length(spike_detection_ds)
    dspike_vector(x,1) = length(spike_detection_ds{x});
end
dspike_vector_max = max(dspike_vector);
    for j = 1:length(spike_detection_ds)
        if length(spike_detection_ds{j}) ~= dspike_vector_max
            dif = abs(length(spike_detection_ds{j}) - dspike_vector_max);
            spike_detection_ds{j}(end+dif) = 0;
        end
        for n = 1:dspike_vector_max
            if spike_detection_ds{j}(n) == 0
                spike_detection_ds{j}(n) = NaN;
            end
        end
        plot(spike_detection_ds{j},j,'.','MarkerSize',10,'Color','black'); hold on;
        dspike_plot(:,j) = cell2mat(spike_detection_ds(j)'); 
        
    end
    axis([0 T(end) 0 length(spike_detection_ds)+1]);
% figure('Name','VmLowpass','NumberTitle','off')
% plot(T,lowpass_Vm); hold on;
% ax = gca;
% ax.Box = 'off';
% figure('Name','dLowpass','NumberTitle','off')
% plot(T(end),ds_lowpass);hold on;
% ax = gca;
% ax.Box = 'off';
figure('Name','Raster Plot','NumberTitle','off')
 for v = 1:length(spike_detection_ds)
     subplot(2,1,2),plot_raster(spike_detection_ds{v},1,0.2,'k',0.25); hold on;
 end
 xlim([0 T(end)]);
 ylim([0 length(spike_detection_ds)+2]);
 ax = gca;
 ax.Box = 'off';
 rastertime = toc
%% Spike Width

% for w = 1:length(spike_time_begin)
%     if length(spike_time_begin{1,w}) ~= 0
%         Spike_width(w) = abs(spike_time_begin{1,w}(:,1)- spike_time_begin{1,w}(:,length(spike_time_begin{1,w})))*1000; %Spiking width in ms
%     end
% end

%%

% 
% %Resting Vm Changes
% figure('Name','Resting Membrane (V)');
% plot(linspace(0,length(V_sub),length(V_sub)),V_sub);
% % Spike Number
% for j = 1:length(spike_time)
%    boxplot((spike_time{1,j})); 
% end
% 
% %%
% figure()
% plot(T, Vm(:,2), 'Linewidth',1); axis off;
% axis([6500 8729 -80 60])
% 
%%