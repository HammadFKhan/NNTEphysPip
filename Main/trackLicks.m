function licks =  trackLicks(lickData)

I = find(lickData==1);
I_diff = diff(I);
idx = find(I_diff~=1);
idx = idx+1;
licks = horzcat(I(1),I(idx));
%VR.licks = licks;