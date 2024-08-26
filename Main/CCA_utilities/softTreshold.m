%%%% Soft tresholding of vector w to make its L1 < c
function [ws,d]=softTreshold(w,c)

for  d=linspace(0,max(w),1000)
    ws=sign(w).* max(abs(w)-d,0);
    ws=ws/norm(ws);
    if sum(abs(ws))<=c
        break;
    end
end
end
