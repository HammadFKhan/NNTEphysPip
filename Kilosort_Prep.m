%% Kilosort data preprocessing

clear; clc; close all;
read_Intan_RHD2000_file

datI = int16(amplifier_data);

fid = fopen(strcat('KilosortData/',filename(1:end-4)),'w');
fwrite(fid,datI,'int16');
fclose(fid);
