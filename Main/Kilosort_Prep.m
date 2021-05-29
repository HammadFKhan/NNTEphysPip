%% Kilosort data preprocessing

clear; clc; close all;

datI = int16(amplifier_data);

fid = fopen(strcat('KilosortData/',filename(1:end-4)),'w');
fwrite(fid,datI,'int16');
fclose(fid);
