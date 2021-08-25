function filtData = fastpreprocess_filtering(IntanData,Fs)

if nargin < 3 || strcmp(Fs,'')
    Fs = 20000;
    disp(['Sampling rate set at ' num2str(Fs) ' Hz for filtering']);
end
% Notch Filtering
channel_num = size(IntanData,1);
d = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    'DesignMethod','butter','SampleRate',Fs);

d1 = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',119,'HalfPowerFrequency2',121, ...
    'DesignMethod','butter','SampleRate',Fs);
H = waitbar(0,'Filtering...');
for i = 1:channel_num
waitbar(i/channel_num(end),H)
rawData60 = filtfilt(d,double(IntanData(i,:)));
rawData120 = filtfilt(d1,rawData60);  
filtData.rawData(i,:) = rawData120;
filtData.channel_num = channel_num;
end
close(H)

