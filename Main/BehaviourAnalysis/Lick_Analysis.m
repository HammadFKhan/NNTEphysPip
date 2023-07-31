%Licks per Trial
addpath(genpath('main'));
pathname = uigetdir(pwd,'Input Directory');
pathname = fullfile(pathname);

directory = dir(pathname);
count = 1;
for idx = 3:length(directory)
  file = directory(idx).folder;
  path = directory(idx).name;
  Intan = read_Intan_RHD2000_file(file,path);  
  licks = trackLicks(Intan.board_dig_in_data(3,:));
  triggers = trackLicks(Intan.board_dig_in_data(2,:));
  allLicks(count) = length(licks);
  allTriggers(count) = length(triggers);
  count = count+1;
end % load Intan files

%%
for i = 1:size(VR_data.AvgVel,2)
    avgVelt(i) = mean(VR_data.AvgVel{1,i},1)
end
