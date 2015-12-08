function [ weights ] = generate_weights(binary_masks)

% Find number of pixels (area) in each bin
window_radius = 40;
num_a_bins = 12;
num_r_bins = 5;
num_pts = size(binary_masks, 3);
total_area = zeros(num_r_bins, num_a_bins);

% Preliminary calculations

num_ones = zeros(num_r_bins, num_a_bins, num_pts);
angle_width = 2*pi/num_a_bins;
dist_width = log(window_radius)/num_r_bins;
% Put ones in bin
for r = -window_radius:window_radius
    for c = -window_radius:window_radius
        d_bin = max(1, ceil(1/2*log(r*r+c*c)/dist_width));
        if (d_bin <= num_r_bins)
            a_bin = ceil((atan2(r, c)+pi)/angle_width);
            r_win = r + window_radius+1;
            c_win = c + window_radius+1;
            num_ones(d_bin, a_bin, :) = num_ones(d_bin, a_bin, :) +...
                                        binary_masks(r_win, c_win, :);
            total_area(d_bin, a_bin) = total_area(d_bin, a_bin) + 1;
        end
    end
end
%apply angular blur
Gx = normpdf(-1:1:1, 0, 0.7);
total_area = conv2(padarray(total_area, [0, 1], 'circular'), Gx, 'valid');

weights = num_ones ./ repmat(total_area, 1, 1, num_pts);

end

