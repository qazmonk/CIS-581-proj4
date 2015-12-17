function [ warped_pts ] = apply_tps(a1_x, ax_x, ay_x, w_x, a1_y, ...
                                    ax_y, ay_y, w_y, ctr_pts, pts)

U = pdist2(pts, ctr_pts);
U = spfun(@(r) (r.^2) .* log(r.^2), U);
M = [ones(length(pts), 1), pts(:, 1), pts(:, 2), U];
warped_pts = [M*[a1_x;ax_x;ay_x;w_x], M*[a1_y;ax_y;ay_y;w_y]];
end