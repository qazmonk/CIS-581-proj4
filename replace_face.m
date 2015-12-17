function [ replaced ] = replace_face(I, face, mask, icen, fcen, ipts, fpts)

num_pts = length(ipts);
fpts = fpts - repmat(fcen, num_pts, 1);
ipts = ipts - repmat(icen, num_pts, 1);
[fi_a1_x, fi_ax_x, fi_ay_x, fi_w_x] = est_tps(fpts, ipts(:, 1));
[fi_a1_y, fi_ax_y, fi_ay_y, fi_w_y] = est_tps(fpts, ipts(:, 2));
[if_a1_x, if_ax_x, if_ay_x, if_w_x] = est_tps(ipts, fpts(:, 1));
[if_a1_y, if_ax_y, if_ay_y, if_w_y] = est_tps(ipts, fpts(:, 2));



[face_nr, face_nc, ~] = size(face);
[face_xs, face_ys] = meshgrid(1:face_nc, 1:face_nr);
face_inds = sub2ind([face_nr, face_nc], face_ys(:), face_xs(:));
hull = bwmorph(mask, 'remove');
face_inds(hull(face_inds) == 0) = [];
[face_ys, face_xs] = ind2sub([face_nr, face_nc], face_inds);

num_hull_pts = length(face_xs);
face_xs = face_xs - repmat(fcen(1), num_hull_pts, 1);
face_ys = face_ys - repmat(fcen(2), num_hull_pts, 1);


warped_hull = apply_tps(fi_a1_x, fi_ax_x, fi_ay_x, fi_w_x, fi_a1_y, ...
                        fi_ax_y, fi_ay_y, fi_w_y, fpts, [face_xs ...
                    face_ys]);
warped_hull = warped_hull + repmat(icen, num_hull_pts, 1);

K = convhull(warped_hull);

[I_nr, I_nc, ~] = size(I);
[I_xs, I_ys] = meshgrid(1:I_nc, 1:I_nr);

in_hull = inpolygon(I_xs(:), I_ys(:), warped_hull(K, 1), warped_hull(K, 2));
I_xs(~in_hull) = [];
I_ys(~in_hull) = [];

num_I_pts = length(I_xs);
I_xs = I_xs(:) - repmat(icen(1), num_I_pts, 1);
I_ys = I_ys(:) - repmat(icen(2), num_I_pts, 1);

warped_I_pts = apply_tps(if_a1_x, if_ax_x, if_ay_x, if_w_x, if_a1_y, ...
                         if_ax_y, if_ay_y, if_w_y, ipts,...
                         [I_xs,I_ys]); 
warped_I_pts = warped_I_pts + repmat(fcen, num_I_pts, 1);

num_inds = length(I_xs);
ch_inds = [repmat(1, num_inds, 1);
           repmat(2, num_inds, 1);
           repmat(3, num_inds, 1);];

I_xs = repmat(round(I_xs + repmat(icen(1), num_I_pts, 1)), 3, 1);
I_ys = repmat(round(I_ys + repmat(icen(2), num_I_pts, 1)), 3, 1);
I_inds = sub2ind([I_nr I_nc 3], I_ys(:), I_xs(:), ch_inds);

wI_xs = repmat(round(warped_I_pts(:,1)), 3, 1);
wI_ys = repmat(round(warped_I_pts(:,2)), 3, 1);
wI_inds = sub2ind([face_nr face_nc 3], wI_ys(:), wI_xs(:), ch_inds);



replaced = I;
replaced(I_inds) = face(wI_inds);

end