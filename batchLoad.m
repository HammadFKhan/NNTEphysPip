%% Parallel Batch processing 
function [pathname,VR_data] = batchLoad()
pathname = strcat(uigetdir(pwd,'Input Directory'),'\');
directory = dir(pathname);
for idx = 1:length(directory)
    fullfile = strcat(directory(idx).folder,'\',directory(idx).name)
    try
        VR_data{idx} = table2array(readtable(fullfile));
    catch ME
        VR_data{idx} = NaN;
        continue
    end
    
    
    
end

end

