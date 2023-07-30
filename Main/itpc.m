%% Compute and plot TF-ITPC for one electrode
function tf = itpc(lfp,timestamps,Fs)
% add argument for specific electrode!!
% if nargin<3
%     numelectrode = [];
% else
%     numelectrode = 64;
% end
% wavelet parameters
num_frex = 30; % Frequency resolution
min_freq =  5; % Lower frequency bound
max_freq = 40; % Upper frequency bound
timestamps = timestamps.*Fs;
pnts = timestamps(1,2)-timestamps(1,1); % Sample size of each trial
trials = length(timestamps); % Running event trials

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
%%
% build data format time x trials
for i = 1:trials
    try
        data(:,i) = lfp(ceil(timestamps(i,1)):ceil(timestamps(i,2)));
    catch ME
        continue;
    end
end

% FFT of data (doesn't change on frequency iteration)
dataX = fft( reshape(data,[],1) ,nConv)';

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
figure(), clf
contourf(1:pnts,frex,tf,10,'linecolor','none')
set(gca,'clim',[0 .5],'ydir','normal')
title('ITPC')
colormap(jet),colorbar
