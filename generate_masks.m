function [ masks new_locs] = generate_masks( mask, locs )
%GENERATE_MASKS Summary of this function goes here
%   Detailed explanation goes here


s = regionprops(mask, 'centroid');
center = s.Centroid;
new_locs = locs - repmat(center, length(locs), 1);
masks = get_windows(mask, locs, 40);

end

