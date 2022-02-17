% Concatenate trials

addpath(genpath('main'));
pathname = uigetdir(pwd,'Input Directory');
pathname = fullfile(pathname);
directory = dir(fullfile(pathname,'*.rhd')); %Parses RHD files
count = 1;
downsampleRate = 4;
targetedFs = 8192;
L = length(directory);
for idx = 1:L
    file = directory(idx).folder;
    path = directory(idx).name;
    Intan = read_Intan_RHD2000_file(file,path); 
    Fs =  Intan.frequency_parameters.amplifier_sample_rate;
    allIntan{count} = resample(Intan.amplifier_data',targetedFs,Fs);
    count = count+1;
end % load Intan files
Fs = Intan.frequency_parameters.amplifier_sample_rate;
% Concatenate intan files for the whole session
disp('Combining...')
Intan.allIntan = vertcat(allIntan{:})';
disp('Compressing...')
Intan.allIntan = single(Intan.allIntan);
% Adjust electrode order by depth
UCLA_probe_map
Intan.allIntan  = Intan.allIntan(s.sorted_electrodes,:);
% Fix recording offset
Intan.offset = 1; % second
Intan.offsetSample = Intan.frequency_parameters.amplifier_sample_rate*Intan.offset;
disp(['Adjusting for ' num2str(Intan.offset) ' second offset']);
Intan.allIntan = Intan.allIntan(:,Intan.offsetSample:size(Intan.allIntan,2));
clear amplifier_data t_amplifier frequncy_parameters notes aux_input_channels...
    aux_input_data board_dig_in_channels board_dig_in_data amplifier_channels


