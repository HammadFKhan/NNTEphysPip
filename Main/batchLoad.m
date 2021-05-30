%% Parallel Batch processing 
function [pathname,VR_data] = batchLoad()
pathname = uigetdir(pwd,'Input Directory');
pathname = fullfile(pathname);
directory = dir(pathname);
for idx = 1:length(directory)
    filename = fullfile(directory(idx).folder,directory(idx).name)
    try
        VR_data{idx} = table2array(readtable(filename));
    catch ME
        VR_data{idx} = NaN;
        continue
    end
    
    
    
end

end

