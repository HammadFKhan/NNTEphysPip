%% Intan Concatenation/Preprocessing Dual Shanks
% Searches for all active analog/digital inputs for concatenation with
% target resampling. Target resampling will happen at the end of the data
% set for only digital analog I/O pins. To preserve spiking fidelity Intan
% amplifier data is recommended not to resampled and delete channel data
% after LFP has been calculated
% output should be a downsampled amplifer data to 5000 Hz, deleted Intan
% amplifier data and kilosort prepped data file. Then data is send to memory mapped file. This allows for proper
% memory management
% Function call for dual shanks
function ds_filename = intanPreprocessingDualShanks(pathname,chanMapFile,intandsFlag,activeElectrodes)
disp(['Using ' chanMapFile ' as electrode map'])
pause(1)
load(chanMapFile)
pathname = fullfile(pathname);
directory = dir(fullfile(pathname,'*.rhd')); %Parses RHD files
targetedFs = 2000;
L = length(directory);
% Now we build the memory map file if file doesnt exist % only for LFP/behaviour
% data. Spikes are sent to .bin files for kilosort

%because of kilosort we have to make seperate .bin folders for correct execution of data
if ~exist(fullfile(pathname,chanMapFile(1:end-4)),'dir')
    dirFlag = mkdir(pathname,chanMapFile(1:end-4));
    if ~dirFlag
        error('Failed directory creation')
    end
    disp('New directory made for .bin file!')
end
npathname = fullfile(pathname,chanMapFile(1:end-4));
ds_filename = fullfile(npathname,['intan_ds_data_',chanMapFile(1:end-4),'.mat']); % check incremented file name for recordings
kilosort_filename = fullfile(npathname,['kilosort_',chanMapFile(1:end-4),'.bin']);
if exist(ds_filename,'file') %check if downsampled data file already exists
    warning('Preprocessed file already exists! Data will now be overrided')
end

data = matfile(ds_filename,'Writable',true);
kname = ['kilosort_',chanMapFile(1:end-4)];
% Grab parameters and generate structures based on first file
idx = 1;
path = directory(idx).folder;
file = directory(idx).name;
Intan = read_Intan_RHD2000_file(path,file);
data.Fs =  Intan.frequency_parameters.amplifier_sample_rate;
data.path = path;
intanOffset = 1;
% Removing first second of data
disp(['Adjusting for ' num2str(intanOffset) ' second offset']);
Intan.amplifier_data = Intan.amplifier_data(activeElectrodes,data.Fs*intanOffset:size(Intan.amplifier_data,2));
Intan.t_amplifier = Intan.t_amplifier(:,data.Fs*intanOffset:size(Intan.t_amplifier,2));

% running kilosort prep file
if ~exist(kilosort_filename,'file')
    temp = Intan.amplifier_data(s.sorted_electrodes,:); % sort electrodes since we use sorted electrodes in kilosort
    kilosortPrep2(temp,path,kname)
else
    warning('An existing kilosort.bin file exists! Deleting existing kilosort version')
    pause(1)
    delete(kilosort_filename)
    temp = Intan.amplifier_data(s.sorted_electrodes,:);
    kilosortPrep2(temp,path,kname)
end
% Now downsample data for LFP
amplifierData{idx} = resample(Intan.amplifier_data',targetedFs,data.Fs)';
amplifierTime{idx} = downsample(Intan.t_amplifier',round(data.Fs/targetedFs),1)';

if ~isempty(Intan.board_dig_in_data) % Checks for digital traces
    Intan.board_dig_in_data = Intan.board_dig_in_data(:,data.Fs*intanOffset:size(Intan.board_dig_in_data,2));
    digitalChannels{idx} = downsample(Intan.board_dig_in_data',round(data.Fs/targetedFs),1)';
    data.digitalChannelsinfo = Intan.board_dig_in_channels; % save meta data (do once)
end
if ~isempty(Intan.board_adc_data) % Checks for analog traces
    Intan.board_adc_data = Intan.board_adc_data(:,data.Fs*intanOffset:size(Intan.board_adc_data,2));
    analogChannels{idx} = downsample(Intan.board_adc_data',round(data.Fs/targetedFs),1)';
    data.analogChannelsinfo = Intan.board_adc_channels; % save meta data (do once)
end
for idx = 2:L
    path = directory(idx).folder;
    file = directory(idx).name;
    Intan = read_Intan_RHD2000_file(path,file);
    Intan.amplifier_data = Intan.amplifier_data(activeElectrodes,:);
    if idx== L %subtract the last second off the recording
        disp(['Adjusting for ' num2str(intanOffset) ' second offset']);
        Intan.amplifier_data = Intan.amplifier_data(:,1:(size(Intan.amplifier_data,2)-data.Fs*intanOffset));
        Intan.t_amplifier = Intan.t_amplifier(:,1:(size(Intan.t_amplifier,2)-data.Fs*intanOffset));
        if exist('digitalChannels','var')
            Intan.board_dig_in_data = Intan.board_dig_in_data(:,1:(size(Intan.board_dig_in_data,2)-data.Fs*intanOffset));
        end
        if exist('analogChannels','var')
            Intan.board_adc_data = Intan.board_adc_data(:,1:(size(Intan.board_adc_data,2)-data.Fs*intanOffset));
        end
    end
    temp = Intan.amplifier_data(s.sorted_electrodes,:);
    kilosortPrep2(temp,path,kname)
    amplifierData{idx} = resample(Intan.amplifier_data',targetedFs,data.Fs)';
    amplifierTime{idx} = downsample(Intan.t_amplifier',round(data.Fs/targetedFs),1)';
    if exist('digitalChannels','var')
        digitalChannels{idx} = downsample(Intan.board_dig_in_data',round(data.Fs/targetedFs),1)';
    end
    if exist('analogChannels','var')
        analogChannels{idx} = downsample(Intan.board_adc_data',round(data.Fs/targetedFs),1)';
    end
        
end
% Combine cells and save into data
if intandsFlag % If we want to make an intan file (save time for hdd loading
    fprintf('Compressing data...')
    amplifierData = horzcat(amplifierData{:});
    fprintf('done\n')
    % sort electrodes
    fprintf('Saving amplifier data...')
    amplifierData = amplifierData(s.sorted_electrodes,:);
    data.amplifierData = amplifierData;
    fprintf('Saving everything else...')
    data.chanMapFile = chanMapFile;
    if exist('digitalChannels','var')
        data.digitalChannels = horzcat(digitalChannels{:});
    end
    if exist('analogChannels','var')
        data.analogChannels = horzcat(analogChannels{:});
    end
    data.amplifierTime = horzcat(amplifierTime{:});
    fprintf('done\n')
end
data.targetedFs = targetedFs;
fpath = fileparts(ds_filename);
data.fpath = fpath;
clearvars -except ds_filename
end