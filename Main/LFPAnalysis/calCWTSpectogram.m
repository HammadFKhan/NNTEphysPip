function [powerCWT, fwt] = calCWTSpectogram(x,t,Fs,VoicesperOctave,flimit,showplot,baselinesub)

%Ref - https://www.cosmos.esa.int/documents/1655127/1655136/Torrence_1998_Wavelet_Guide_BAMS.pdf/001d8327-b255-3024-a2f0-ce02e33ac98f
if ~exist('baselinesub','var')
    BLsub = 0;
else
    BLsub = baselinesub;
    nBLSub = 0.05*Fs; % in number of points
end

assert( numel(x)==numel(t), 'Number of elements in x and t dont match' );

[wt,fwt] = cwt(x,'amor',Fs,'VoicesperOctave',VoicesperOctave,'FrequencyLimits',flimit); % Calculate the CWT
powerCWT = (abs(wt).^2)/abs(var(x,1)); % Estimate power from CWT 
% The power calculated is normalzied/relative to th ewhite noise power 

if BLsub == 1
    BLPower = mean(powerCWT,2);
    powerCWTBLsub = bsxfun(@rdivide,powerCWT,BLPower);
    powerCWT = powerCWTBLsub;
end

if showplot
    plotSpectrogram(10*log10(powerCWT),t,fwt,'contourf','Wavelet Based Spectrogram','Time (s)','Frequency (Hz)')
end


