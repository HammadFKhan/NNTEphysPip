function filtData = preprocess_filtering(IntanData,time,Fs)
if nargin < 3 || strcmp(Fs,'')
    Fs = 20000;
    disp(['Sampling rate set at ' num2str(Fs) ' Hz for filtering']);
end
channel_num = size(IntanData,1);
d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',Fs);
% Design butterworth filters for single unit
Fc = [300 3000];
Wn = Fc./(Fs/2);
[b1,a1] = butter(6,Wn,'bandpass');

% Design butterworth filters for LFP
Fc = [25];
Wn = Fc./(Fs/2);
[b4,a4] = butter(2,Wn,'low');

% Design butterworth filter for Ripples
Fc = [140 200];
Wn = Fc./(Fs/2);
[b3,a3]=butter(4,Wn,'bandpass');

H = waitbar(0,'Filtering...');
for i = 1:channel_num
waitbar(i/channel_num(end),H)
rawData = filtfilt(d,double(IntanData(i,:)));
lowpassData  = filtfilt(b4,a4,rawData);
bandpassData = filtfilt(b1,a1,rawData);

pow = fastrms(lowpassData,15);    
    mRipple(i) = mean(pow);
    meRipple(i) = median(pow);
    mmRippleRatio(i) = mRipple(i)./meRipple(i);
    
filtData.rawData(i,:) = rawData;
filtData.lowpassData(i,:) = lowpassData;
filtData.bandpassData(i,:) = bandpassData;
% filtData.LFP.(['Channel' num2str(i)]) = lowpassData;
% filtData.MUA.(['Channel' num2str(i)]) = bandpassData;
filtData.ripple.data(:,i) = FiltFiltM(b3,a3,rawData);
filtData.channel_num = channel_num;
end

medianSignal = median(filtData.bandpassData);
commonModeAvg = filtData.bandpassData-medianSignal;
filtData.commonModeAvg = commonModeAvg;
mmRippleRatio(mRipple<1) = 0;
mmRippleRatio(meRipple<1) = 0;

[~,loc] = max(mmRippleRatio);
bestLFPchan = loc;

disp(['Recommended LFP channel ' bestLFPchan]);
filtData.ripple.timestamps(:,1) = time;
filtData.bestLFPchan = bestLFPchan;
filtData.LFPpower = mmRippleRatio;
close(H)

