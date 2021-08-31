function LFP = fastpreprocess_filtering(intanData,Fs)

if nargin < 2 || strcmp(Fs,'')
    Fs = 20000;
    disp(['Sampling rate set at ' num2str(Fs) ' Hz for filtering']);
end

LFP.Fs = Fs;
downSampleFreq = 1024;

% Notch Filtering
channel_num = size(intanData,1);
d = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    'DesignMethod','butter','SampleRate',downSampleFreq);

d1 = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',119,'HalfPowerFrequency2',121, ...
    'DesignMethod','butter','SampleRate',downSampleFreq);

Fc = [254];
Wn = Fc./(Fs/2);
b = fir1(10,Wn,'low');

% Design butterworth filters for LFP
% Fc = [128];
% Wn = Fc./(Fs/2);
% [b4,a4] = fir1(2,Wn,'low');
disp(['Downsampling to ' num2str(downSampleFreq) ' Hz...'])
downsample_Data = resample(double(intanData)',downSampleFreq,Fs)';
disp('Filtering...')
lowpass = zeros(channel_num,size(downsample_Data,2));
rawData60 = filtfilt(d,downsample_Data);
rawData120 = filtfilt(d1,rawData60); 
LFP.LFP = filtfilt(b,1,rawData120);
LFP.channel_num = channel_num;

