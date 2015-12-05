function [ windows ] = get_windows( I, pts, window_radius )
%GETWINDOWS Summary of this function goes here
%   Detailed explanation goes here

%pad edges and shift locations accordingly
num_pts = size(pts, 1);
I_pad = double(padarray(I, [window_radius, window_radius], 'symmetric'));
shifted_pts = pts + repmat([window_radius, window_radius], num_pts, 1);

window_size = [2*window_radius+1, 2*window_radius+1];
windows = zeros([window_size, num_pts]);
for r = -window_radius:window_radius
    for c = -window_radius:window_radius
        pt = shifted_pts + repmat([r, c], num_pts, 1);
        pt_inds = sub2ind(size(I_pad), pt(:, 1), pt(:, 2));
        r_win = r + window_radius+1;
        c_win = c + window_radius+1;
        windows(r_win, c_win, :) = reshape(I_pad(pt_inds), 1, 1, num_pts);
    end
end


end

