function [ output_args ] = find_face( I, poses )
%FIND_FACE Summary of this function goes here
%   Detailed explanation goes here

%Generate codebook from faces (use function generate_codebook)
%generate descriptors for interest point in each pose image
    % Search grid for points near edge point
    % generate HOG or shape context feature for point along with position
    % and mask information

%for each scale we want to search
    %generate interest points in <I> using grid search again
    %generate HOG or shape context features and find k closest matches from
    %codebook. Use matches to create score map
%do non maximal supresssion across scales to find hypothesis points

%Pick best hypothesis and crop out face and replace it with pose that
%generated the most votes. Come up with some correspondance between points
%then warp and blend
end

