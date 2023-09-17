%% Intan Concatenation/Preprocessing V2
% Searches for all active analog/digital inputs for concatenation with
% target resampling. Target resampling will happen at the end of the data
% set for only digital analog I/O pins. To preserve spiking fidelity Intan
% amplifier data is recommended not to resampled and delete channel data
% after LFP has been calculated
% output should be a downsampled amplifer data to 5000 Hz, deleted Intan
% amplifier data and kilosort prepped data file. Then data is send to memory mapped file. This allows for proper
% memory management
function ds_filename = intanPreprocessing2
addpath(genpath('Main'));
chanMapFile = 'UCLA_chanmap_fixed.mat';
disp(['Using ' chanMapFile ' as electrode map'])
pause(1)
load(chanMapFile)
pathname = uigetdir(pwd,'Input Directory');
pathname = fullfile(pathname);
directory = dir(fullfile(pathname,'*.rhd')); %Parses RHD files
targetedFs = 2000;
L = length(directory);
% Now we build the memory map file if file doesnt exist % only for LFP/behaviour
% data. Spikes are sent to .bin files for kilosort
ds_filename = fullfile(pathname,'intan_ds_data.mat'); % check incremented file name for image
kilosort_filename = fullfile(pathname,'kilosort.bin');
if exist(ds_filename,'file') %check if downsampled data file already exists
    error('Preprocessed file already exists! Please remove from directory')
end

data = matfile(ds_filename,'Writable',true);
% Grab parameters and generate structures based on first file
idx = 1;
path = directory(idx).folder;
file = directory(idx).name;
Intan = read_Intan_RHD2000_file(path,file);
data.Fs =  Intan.frequency_parameters.amplifier_sample_rate;
data.path = path;
intanOffset = 1;
% Removing first second of data6
disp(['Adjusting for ' num2str(intanOffset) ' second offset']);
Intan.amplifier_data = Intan.amplifier_data(:,data.Fs*intanOffset:size(Intan.amplifier_data,2));
Intan.t_amplifier = Intan.t_amplifier(:,data.Fs*intanOffset:size(Intan.t_amplifier,2));

% running kilosort prep file
if ~exist(kilosort_filename,'file')
    kilosortPrep2(Intan.amplifier_data,path)
else
    warning('An existing kilosort.bin file exists! Deleting existing kilosort version')
    pause(1)
    delete(kilosort_filename)
    kilosortPrep2(Intan.amplifier_data,path)
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
    if idx==L %subtract the last second off the recording
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
    kilosortPrep2(Intan.amplifier_data,path)
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
fprintf('Compressing data...')
amplifierData = horzcat(amplifierData{:});
fprintf('done\n')
% sort electrodes
fprintf('Saving amplifier data...')
amplifierData = amplifierData(s.sorted_electrodes,:); 
data.amplifierData = amplifierData;
fprintf('Saving everything else...')
data.chanMapFile = chanMapFile;
data.digitalChannels = horzcat(digitalChannels{:});
data.analogChannels = horzcat(analogChannels{:});
data.amplifierTime = horzcat(amplifierTime{:});
fprintf('done\n')
data.targetedFs = targetedFs;
data.fpath = Intan.path;
clearvars -except ds_filename
end

