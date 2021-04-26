%% Kilosort data preprocessing

clear; clc; close all;
read_Intan_RHD2000_file

datI = int16(amplifier_data(5:36,:));

fid = fopen('myNewFile.bin','w');
fwrite(fid,datI,'int16');
fclose(fid);
