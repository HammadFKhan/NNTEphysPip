function output = preprocessSpike(data,Fs)
Fc = [300 2000];
Wn = Fc./(Fs/2);
b = fir1(5000,Wn,'bandpass');
rawspikeTrace = filtfilt(b,1,double(data)');
rawspikeTrace = rawspikeTrace';
commonModeAvg = rawspikeTrace-mean(rawspikeTrace);
whitenedSpikeTrace = commonModeAvg;
output.rawspikeTrace = rawspikeTrace;
output.whitenedSpikeTrace = whitenedSpikeTrace;
end