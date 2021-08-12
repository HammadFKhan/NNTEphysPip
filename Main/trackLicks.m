function licks =  trackLicks(lickData)

I = find(lickData==1);
I_diff = diff(I);
idx = find(I_diff~=1);
idx = idx+1;
try
    licks = horzcat(I(1),I(idx));
catch ME
    licks = 0;
end

%VR.licks = licks;