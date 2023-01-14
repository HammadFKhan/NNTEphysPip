function [CSD] = CustomCSD(LFP,spacing,rate,plotCSD,smoothCSD)
%% a custom CSD gathered from multiple sources for a clean CSD plot
% Inputs:
%   LFP = LFP in uV. Rows are recording sites, columns are time
%   spacing = electrode z-sspacing in microns
%   rate = sampling rate
%   plotCSD = true/false. true generates a new CSD figure
%   smoothCSD = true/false. true will interpolate the CSD in the Z
%               direction by a factor of 3
 
%convert to volts and meters
spacing = spacing/1e6;
data    =  LFP./1e6;

%% 4 add vaknin electrodes
%add the top electrode two more time and bottom electrode two time for a
%total of 4 Vaknin electrodes. Final CSD will equal the number of real
%electrodes
vakData = zeros(size(data,1)+4,size(data,2));
vakData(1,:) = data(1,:);
vakData(2,:) = data(1,:);
vakData(3:3+size(data,1)-1,:) = data;
vakData(end-1,:) = data(end,:);
vakData(end,:) = data(end,:);

%% hamming filter
%currently not necessary with n = 2 in CSD second derivative approx
 
%3 pt hamming filter
% for zz = 2:size(vakData,1)-1
%     hamFiltData(zz,:) = 0.23.*vakData(zz+1,:) + 0.23.*vakData(zz-1,:) + 0.54.*vakData(zz,:);
% end
% vakData = hamFiltData;

%% compute CSD
%approx of second spatial derivative in Freemon and Nicholson
% i_m = (V_b + V_a - 2*V_o)/(2*delta_z)*2
% where
% V_o is the center electrode
% V_b is the upper electrode
% and V_a is the lower electrode
% "upper and lower" can differ depeding on the equation, but it is less
% noisy if the center electrode is z = 0 and upper/lower at z +/- 2.
% Thus, n = 2.
n = 2;

%allocate CSD memory
CSD = repmat(NaN,size(vakData,1),size(vakData,2));

%do the calculation from 3rd electrode to 3rd to last electrode
for ii = 3:size(vakData,1)-2 % 2nd to 2nd last electode
    %top electrode potential
    V_b = vakData(ii+1*n,:);
    %center electrode potential
    V_o = vakData(ii,:);
    %bottom electrode potential
    V_a = vakData(ii-1*n,:);
    
    %CSD
    CSD(ii,:) = (V_b+V_a-2*V_o)./(2*spacing)^2;
end

%now remove the vaknin electrodes
CSD(size(vakData,1)-1:size(vakData,1),:) = [];
CSD(1:2,:) = [];
% CSD = -1.*CSD;

if(smoothCSD)
    for ii = 1:size(CSD,2)
        CSD_Smooth(:,ii) = interp(CSD(:,ii),5);
    end   
    CSD  = CSD_Smooth;
end

if(plotCSD)
    t = (1:size(CSD,2))./rate.*1000;
    depth = ([1:size(CSD,1)]-1)*spacing/5*1e6;
    figure
    imagesc(t,depth,-1.*CSD)
    colormap(jet)
%     caxis([-1E4,1E4])
end
