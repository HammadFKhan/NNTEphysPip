function intanData = preprocess_filtering(IntanData)
channel_num = size(IntanData,1);
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
H = waitbar(0,'Filtering...');
for i = 1:channel_num
waitbar(i/channel_num(end),H)
rawData = filtfilt(d,IntanData(i,:));
lowpassData  = filtfilt(b4,a4,rawData);
bandpassData = filtfilt(b1,a1,rawData);

intanData.rawData(i,:) = rawData;
intanData.lowpassData(i,:) = lowpassData;
intanData.bandpassData(i,:) = bandpassData;

intanData.LFP.(['Channel' num2str(i)]) = lowpassData;
intanData.MUA.(['Channel' num2str(i)]) = bandpassData;

intanData.channel_num = channel_num;
end
close(H)

