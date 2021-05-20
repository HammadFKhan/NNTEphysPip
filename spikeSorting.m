function Spikes = spikeSorting(intanData,Fs)
if nargin < 2 || strcmp(Fs,'')
    Fs = 20000;
    disp(['Sampling rate set at ' num2str(Fs) ' Hz for spike detection']);
end
findSpikes = true;
threshStd = 5;
channel_num = intanData.channel_num;
bandpassData = intanData.bandpassData;
spikeCount = zeros(channel_num,1);
H = waitbar(0,'Spike Sorting...');
Spikes.binary = Spike_Detector_Single(bandpassData,2.5,0.3);
for i = 1:channel_num
    waitbar(i/channel_num(end),H)
    thresh = threshStd*std(bandpassData(i,:));
    [pks,locs] = findpeaks(-bandpassData(i,:),Fs,'MinPeakHeight',thresh,'MinPeakDistance',.001);
    %save the spike times for this channel
    window = .00175*Fs;%35ms
    if(~isempty(pks))
        Spikes.allLocs.(['Channel' num2str(i)]) = locs;
        for pp = 1:length(locs)
            thisLoc = locs(pp)*Fs;
            if(thisLoc-window>1 && thisLoc+window<length(bandpassData(i,:)))%%peak can't occur at the very end or very beginning of the data set
                spikeCount(i) = spikeCount(i)+1;
                %extract a 4ms window around the spike peak
                Spikes.allSpikes.(['Channel' num2str(i)])(spikeCount(i),:) = bandpassData(i,thisLoc-window:thisLoc+window);
            end
        end
    else
        Spikes.allLocs.(['Channel' num2str(i)]) = [];
    end
    Spikes.spikeCount = spikeCount;
end
close(H)