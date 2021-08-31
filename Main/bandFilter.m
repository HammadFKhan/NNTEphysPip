function LFP = bandFilter(LFP)

Fs = LFP.Fs;
% Filters
% Design butterworth filters for LFP
Fc = [1 10];
Wn = Fc./(Fs/2);
[b1,a1] = butter(2,Wn,'bandpass');

Fc = [10 30];
Wn = Fc./(Fs/2);
[b2,a2] = butter(2,Wn,'bandpass');

Fc = [30 80];
Wn = Fc./(Fs/2);
[b3,a3] = butter(2,Wn,'bandpass');
disp('Band Filtering...')
LFP.theta_band = filtfilt(b1,a1,LFP.LFP);
LFP.beta_band = filtfilt(b2,a2,LFP.LFP);
LFP.gamma_band = filtfilt(b3,a3,LFP.LFP);
