function [ masks ] = generate_masks( mask, locs )
%GENERATE_MASKS Summary of this function goes here
%   Detailed explanation goes here

masks = get_windows(mask, locs, 40);
end

