% Concatenate trials

addpath(genpath('main'));
pathname = uigetdir(pwd,'Input Directory');
pathname = fullfile(pathname);
directory = dir(pathname);

count = 1;
for idx = 3:length(directory)
    file = directory(idx).folder;
    path = directory(idx).name;
    Intan = read_Intan_RHD2000_file(file,path);   
    allIntan{count} = Intan.board_dig_in_data;
    count = count+1;
end % load Intan files
% Concatenate intan files for the whole session
Intan.allIntan = horzcat(allIntan{:});
clear amplifier_data t_amplifier frequncy_parameters notes aux_input_channels...
    aux_input_data board_dig_in_channels board_dig_in_data amplifier_channels


