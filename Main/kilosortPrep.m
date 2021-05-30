function kilosortPrep(amplifier_data,path,filename)
% Kilosort data preprocessing
datI = int16(amplifier_data);
kilosortOut = path;
[SUCCESS,~,~] = mkdir(kilosortOut,filename(1:end-4));
newDIR = dir(strcat(kilosortOut,filename(1:end-4)));
newDIR = newDIR.folder;
fid = fopen(strcat(newDIR,'/',filename(1:end-4),'.bin'),'w');
fwrite(fid,datI,'int16');
fclose(fid);
disp('Kilosort Conversion Successful!')