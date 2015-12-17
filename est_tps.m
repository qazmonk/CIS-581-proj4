function [a1,ax,ay,w] = est_tps(ctr_pts, target_value)
    p = length(ctr_pts);
    K = pdist2(ctr_pts, ctr_pts);
    K = spfun(@(r) (r.^2).*log(r.^2), K);
    P = [ctr_pts, ones(p, 1)];
    O = zeros(3,3);
    M = [K, P; P', O];
    res = (M+eps*eye(p+3,p+3))\[target_value;0;0;0];
    w = res(1:p);
    ax = res(p+1);
    ay = res(p+2);
    a1 = res(p+3);
end