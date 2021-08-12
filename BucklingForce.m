% Buckling Force Model

t = [6*10^-6,10*10^-6:5*10^-6:25*10^-6]; %Thickness of the shank 16 microns

figure
for i = 1:length(t)
    w = 250*10^-6; %Width of the shank 170 microns
    L = 1*10^-3:0.000025:3*10^-3; % Length in millimeter
    E = 3.2*10^9; %Young's modulus of the shank
    
    top  = pi*pi*E*w*t(i)^3;
    
    bottom = 5.88.*L.^2;
    
    ForceB = top./bottom;
    plot(L*1000,ForceB*1000,'LineWidth',2);hold on
end
title('Buckling Force with Parylene Only')
ylim([0 10]);
yline(1,'r--','Critical Point');
legend('6um','10um','15um','20um','25um')
xlabel('Shank Length (mm)')
ylabel('Force mN')
%% Current Design with just parylene
figure
t = [6*10^-6]; %Thickness of the shank 16 microns
w = 250*10^-6; %Width of the shank 170 microns
L = 1*10^-3:0.000025:3*10^-3; % Length in millimeter
E = 3.2*10^9; %Young's modulus of the shank
top  = pi*pi*E*w*t^3;

bottom = 5.88.*L.^2;

ForceB = top./bottom;
plot(L*1000,ForceB.*1000,'LineWidth',2);hold on
title('Buckling Force with Parylene Only')
xlabel('Shank Length (mm)')
ylabel('Force mN')
ylim([0 1.1]);
yline(1,'r--','Critical Point');
legend('6um')
box off
%%
t = [6*10^-6]; %Thickness of the shank 16 microns
t1 = [5*10^-6:5*10^-6:40*10^-6];
figure
for i = 1:length(t1)
    w = 250*10^-6; %Width of the shank 170 microns
    L = 1*10^-3:0.000025:3*10^-3; % Length in millimeter
    E = 3.2*10^9; %Young's modulus of the shank
    E1 = 4.1*10^9; %Young's Modulus of SU8 
    top  = pi*pi*E*w*t^3;
    top1 = pi*pi*E*w*t1(i)^3;
    bottom = 5.88.*L.^2;
    
    ForceB = top./bottom;
    ForceB1 = top1./bottom;
    total = ForceB+ForceB1;
    hold on
    plot(L*1000,total*1000,'LineWidth',2); hold on;
end
plot(L*1000,ForceB*1000,'LineWidth',2);
title('Buckling Force with Base Parylene and SU8')
ylim([0 30]);
yline(1,'r--','Critical Point');
yline(15,'b--','Dura Penetration');
legend('5um','10um','15um','20um')
xlabel('Shank Length (mm)')
ylabel('Force mN')
