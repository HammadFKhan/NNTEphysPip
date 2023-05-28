function LFP = fastpreprocess_filtering(intanData,Fs)
tic
global useGPU
if nargin < 2 || strcmp(Fs,'')
    Fs = 20000;
    disp(['Sampling rate set at ' num2str(Fs) ' Hz for filtering']);
end

LFP.Fs = Fs;
downSampleFreq = 1024;

% Notch Filtering
channel_num = size(intanData,1);
Fc = 60;
Wn = Fc./(downSampleFreq/2);
QF = Wn/35;
[d,c] = iirnotch(Wn,QF);
Fc = 120;
Wn = Fc./(downSampleFreq/2);
QF = Wn/10;
[d1,c1] = iirnotch(Wn,QF);

Fc = [1 250];
Wn = Fc./(downSampleFreq/2);
b = fir1(5000,Wn,'bandpass');

% Design butterworth filters for LFP
% Fc = [128];
% Wn = Fc./(Fs/2);
% [b4,a4] = fir1(2,Wn,'low');
disp(['Downsampling to ' num2str(downSampleFreq) ' Hz...'])
downsample_Data = resample(double(intanData)',downSampleFreq,Fs)';
disp('Filtering...')
% lowpass = zeros(channel_num,size(downsample_Data,2));


%% GPU case
if useGPU
    buff = gpuArray(downsample_Data);
    datr = filter(d,c,buff); % causal forward filter
    datr = flipud(datr); % reverse time
    datr = filter(d, c, datr); % causal forward filter again
    rawData60 = flipud(datr); % reverse time back
    
    datr1 = filter(d1,c1,rawData60); % causal forward filter
    datr1  = flipud(datr1); % reverse time
    datr1  = filter(d1, c1, datr1); % causal forward filter again
    rawData120 = flipud(datr1); % reverse time back
    
    datr2 = filter(b,1,rawData120); % causal forward filter
    datr2  = flipud(datr2); % reverse time
    datr2  = filter(b, 1, datr2); % causal forward filter again
    lfp = flipud(datr2); % reverse time back
    
    %Output GPU data
    LFP.LFP = gather(lfp);
else
%     rawData60 = filtfilt(d,c,downsample_Data');
%     rawData120 = filtfilt(d1,c,rawData60);
    LFP.LFP = filtfilt(b,1,downsample_Data');
end
LFP.LFP = LFP.LFP';
LFP.channel_num = channel_num;
LFP.downSampleFreq = downSampleFreq;
toc
