function kilosortPrep2(amplifier_data,path,kname)
% Kilosort data preprocessing, appends data now instead of overwriting
if ~exist('kname','var') || isempty(kname)
    kname = 'kilosort';
end
datI = int16(amplifier_data);
kilosortOut = [path,'\',kname(10:end)];
% [SUCCESS,~,~] = mkdir(kilosortOut,filename(1:end-4));
% DIR = dir(strcat(kilosortOut,filename(1:end-4)));
% DIR = DIR.folder;
fid = fopen(strcat(kilosortOut,'\',kname,'.bin'),'a');
fwrite(fid,datI,'int16');
fclose(fid);
disp('Kilosort Conversion Successful!')

