function [waveTCN_behaviorTrace] = getPG(waveTCN_behaviorTrace,parameters)


for ii=1:size(waveTCN_behaviorTrace,2)    
    p = waveTCN_behaviorTrace(ii).xgp;
    wt = waveTCN_behaviorTrace(ii).wt;
    if sum(isnan(p),'all') >= 1
        for ll=1:size(p,3)
            p(:,:,ll) = inpaint_nans(p(:,:,ll),3);
            wt(:,:,ll) = inpaint_nans(wt(:,:,ll),3);
        end
    end
    [pm,pd,dx,dy] = phase_gradient_complex_multiplication(p, parameters.xspacing, parameters.yspacing );
    waveTCN_behaviorTrace(ii).dx = dx;
    waveTCN_behaviorTrace(ii).dy = dy;
    insts = instantaneous_speed(wt,pm);
    [vx, vy] = wavefront_direction(pd,insts);
    waveTCN_behaviorTrace(ii).vx = vx;
    waveTCN_behaviorTrace(ii).vy = vy;
    
end
  
end
