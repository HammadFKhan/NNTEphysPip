%% preprocess all rhd files in a folder
%makes data analysis faster bec i don't have to filter eand threshold for
%spikes every time
%filters data, finds spikes, etc, saves all in mat file

%% setup
clear
%sampling rate
Fs = 20000;

%get names of all intan files
fileList = dir('*.rhd');
numTrials = length(fileList);
% numTrials = 1;
path = [pwd '\'];

%define window around pole movement
windowBefore = 1*Fs;
windowAfter  = 2*Fs;

%load first data set to get number of channels and know sizes of memory
%allocation
read_Intan_RHD2000_fileName([path '\'],fileList(1).name)
numChannels = length(amplifier_channels);
%find beginning of pole movement
poleMoveIndices = find(board_dig_in_data(2,:));
dd = find(diff(board_dig_in_data(2,:)));
touchStartIndex = dd(1);
poleLength = length(dd);
trialLength = length(touchStartIndex-windowBefore:touchStartIndex+windowAfter);

findSpikes = true;
spikeWindow = .00175*Fs;%5ms
threshStd = 5;


%% design filters

%design butterworth filters for single unit
Fc = [250 3000];
Wn = Fc./(Fs/2);
[b1,a1] = butter(6,Wn,'bandpass');

%design butterworth filters for LFP
Fc = [250];
Wn = Fc./(Fs/2);
[b2,a2] = butter(6,Wn,'low');

%% go through each trial, fitler and re-save only the designated window around pole movement

for cc = 1:numChannels
    data.rawData      = zeros(numTrials,trialLength);
    data.lowpassData  = zeros(numTrials,trialLength);
    data.bandpassData = zeros(numTrials,trialLength);
    for tt = 1:numTrials
        %load the first intan file
        read_Intan_RHD2000_fileName(path,fileList(tt).name)

        %find beginning of pole movement
        poleMoveIndices = find(board_dig_in_data(2,:));
        dd = find(diff(board_dig_in_data(2,:)));
        touchStartIndex = dd(1);
        poleLength = length(dd);
        %indices to extract
        trialIndices = touchStartIndex-windowBefore:touchStartIndex+windowAfter;
        
        %get raw data
        rawData = amplifier_data(cc,trialIndices);
        %save raw data in struct
        data.rawData(tt,:) = rawData;

        %get LFP band
        lowpassData  = filtfilt(b2,a2,rawData);
        data.lowpassData(tt,:) = lowpassData;

        %bandpass data for spike band
        bandpassData = filtfilt(b1,a1,rawData);
        data.bandpassData(tt,:) = bandpassData;

        %find spikes
        thresh = threshStd*std(bandpassData);
        [pks,locs] = findpeaks(-bandpassData,Fs,'MinPeakHeight',thresh,'MinPeakDistance',.001);
        %save the spike times for this channel
        if(~isempty(pks))
            data.spikeTimes{tt} = locs;
        else
            data.spikeTimes{tt} = [];
        end
    end
    %save Diginital input 2, which tells us when the pole moves
    data.DigIn2 = board_dig_in_data(2,trialIndices);
    
    save(['Channel' num2str(cc) '-Data.mat'],'data')
end





