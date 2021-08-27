function LFP = downSampleLFP(intanData,Fs)
if nargin < 2 || strcmp(Fs,'')
    disp(['Include sample rate!']);
end
% Low pass filter
Fc = [128];
Wn = Fc./(Fs/2);
b = fir1(2,Wn,'low');

for i = 1:size(intanData,1)
    downsample_LFP = resample(double(intanData(i,:)),1024,Fs);
    filtLFP = filtfilt(b,1,downsample_LFP);
    LFP.downsample(i,:) = filtLFP;
end
