function intanData = spikeSorting(intanData,Fs)
if nargin < 2 || strcmp(Fs,'')
    Fs = 20000;
    disp(['Sampling rate set at' num2str(Fs)]);
end
findSpikes = true;
spikeWindow = .00025*Fs;%5ms
threshStd = 5;
ripples = {};
channel_num = intanData.channel_num;
bandpassData = intanData.bandpassData;
intanData.Spikes.spikeCount = zeros(channel_num,1);
H = waitbar(0,'Filtering...');
for i = 1:channel_num
    waitbar(i/channel_num(end),H)
    thresh = threshStd*std(bandpassData);
    [pks,locs] = findpeaks(-bandpassData,Fs,'MinPeakHeight',thresh,'MinPeakDistance',.001);
    %save the spike times for this channel
    window = .00175*Fs;%35ms
    if(~isempty(pks))
        intanData.Spikes.allLocs.(['Channel' num2str(channel_num(i))]) = locs;
        for pp = 1:length(locs)
            thisLoc = locs(pp)*Fs;
            if(thisLoc-window>1 && thisLoc+window<length(bandpassData))%%peak can't occur at the very end or very beginning of the data set
                intanData.Spikes.spikeCount(channel_num(i)) = spikeCount(channel_num(i))+1;
                %extract a 4ms window around the spike peak
                intanData.Spikes.allSpikes.(['Channel' num2str(channel_num(i))])(spikeCount(channel_num(i)),:) = bandpassData(thisLoc-window:thisLoc+window);
            end
        end
    else
        intanData.Spikes.allLocs.(['Channel' num2str(channel_num(i))]) = [];
    end
end
close(H)