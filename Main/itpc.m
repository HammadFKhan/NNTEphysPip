%% Compute and plot TF-ITPC for one electrode
function tf = itpc(LFP,electrode)
% wavelet parameters
num_frex = 40; % Frequency resolution
min_freq =  2; % Lower frequency bound
max_freq = 40; % Upper frequency bound

pnts = size(LFP.medianLFP,2); % Sample size
trials = size(velocityTrig); % Running event trials

channel2use = electrode; % Set reference channel

% set range for variable number of wavelet cycles
range_cycles = [ 4 10 ];

% other wavelet parameters
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
wavecycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex);
time = -2:1/Fs:2;
half_wave_size = (length(time)-1)/2;

% FFT parameters
nWave = length(time);
nData = pnts*trials; %total size of data (sample points x trials)
nConv = nWave+nData-1;

% FFT of data (doesn't change on frequency iteration)
dataX = fft( reshape(data(channel2use,:,:),1,nData) ,nConv);

% initialize output time-frequency data
tf = zeros(num_frex,pnts);

% loop over frequencies
for fi=1:num_frex
    
    % create wavelet and get its FFT
    s = wavecycles(fi)/(2*pi*frex(fi));
    wavelet  = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));
    waveletX = fft(wavelet,nConv);
    
    % run convolution
    as = ifft(waveletX.*dataX,nConv);
    as = as(half_wave_size+1:end-half_wave_size);
    as = reshape(as,pnts,trials);
    
    % compute ITPC
    tf(fi,:) = abs(mean(exp(1i*angle(as)),2));
end

%% plot results
figure(1), clf
contourf(times,frex,tf,15,'linecolor','none')
set(gca,'clim',[0 .6],'ydir','normal','xlim',[-300 1000])
title('ITPC')
colormap(jet),colorbar