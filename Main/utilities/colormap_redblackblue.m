function ctotal = colormap_redblackblue()
% Colormap to make nice red black blue 
c1 = colormap(flip(hot)); close all;
[grad,~]=colorGradient(c1(end,:),[0 0 1],64);
[grad1,~]=colorGradient([0 0 1],[117 227 255]/256,64);
ctotal = [colormap(flip(hot(128)));grad;grad1];
