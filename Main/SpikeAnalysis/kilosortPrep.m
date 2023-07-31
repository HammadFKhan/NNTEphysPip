function kilosortPrep(amplifier_data,path)
% Kilosort data preprocessing
datI = int16(amplifier_data);
kilosortOut = path;
% [SUCCESS,~,~] = mkdir(kilosortOut,filename(1:end-4));
% DIR = dir(strcat(kilosortOut,filename(1:end-4)));
% DIR = DIR.folder;
fid = fopen(strcat(kilosortOut,'/','kilosort','.bin'),'w');
fwrite(fid,datI,'int16');
fclose(fid);
disp('Kilosort Conversion Successful!')