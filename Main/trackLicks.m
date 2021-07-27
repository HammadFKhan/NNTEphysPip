function VR =  trackLicks(lickData)

lickData  = board_dig_in_data(2,:);
I = find(lickData==1);
I_diff = diff(I);
idx = find(I_diff~=1);
idx = idx+1;
licks = horzcat(I(1),I(idx));
VR.licks = licks;