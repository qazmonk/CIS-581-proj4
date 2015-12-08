function [ sc ] = generate_sc( edges, pts, window_radius)
%GENERATE_SC Summary of this function goes here
%   Detailed explanation goes here

%preliminary calculations for the shape contexts
num_pts = size(pts, 1);
sc = zeros(5, 12, num_pts);
angle_width = 2*pi/12;
dist_width = log(window_radius)/5;

windows = get_windows(edges, pts, window_radius);
%loop over window caluclating distance and angle bins
for r = -window_radius:window_radius
    for c = -window_radius:window_radius
        d_bin = max(1, ceil(1/2*log(r*r+c*c)/dist_width));
        if (d_bin <= 5)
            a_bin = ceil((atan2(r, c)+pi)/angle_width);
            r_win = r + window_radius+1;
            c_win = c + window_radius+1;
            sc(d_bin, a_bin, :) = sc(d_bin, a_bin, :) + ...
                                  windows(r_win, c_win, :);
        end
    end
end

%apply angular blur
Gx = normpdf(-1:1:1, 0, 0.7);
for n = 1:num_pts
    sc(:, :, 1) = conv2(padarray(sc(:, :, n), [0, 1], 'circular'), Gx, 'valid');
end

end

