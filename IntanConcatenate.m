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
    allIntan{count} = Intan.amplifier_data;
    count = count+1;
end % load Intan files
% Concatenate intan files for the whole session
allIntan = horzcat(allIntan{:});


