function [ tps, inlier_ind ] = ransac_est_tps( pts1, pts2, thresh )
%RANSAC_EST_HOMOGRAPHY Summary of this function goes here
%   Detailed explanation goes here

N = 1500;
num_pairs = size(pts1, 1);
votes = cell(N, 1);
max_ind = 1;


for i = 1:N
    sample = datasample(1:num_pairs, 5, 'Replace', false);
    [a1_x,ax_x,ay_x,w_x] = est_tps(pts1(sample, :), pts2(sample, 1));
    [a1_y,ax_y,ay_y,w_y] = est_tps(pts1(sample, :), pts2(sample, 2));

    warped_pts2 = apply_tps(a1_x,ax_x,ay_x,w_x,...
                            a1_y,ax_y,ay_y,w_y, pts1(sample, :), pts2);
    dists = sum((warped_pts2 - pts2).^2, 2);
    %figure('Name', 'Distance Hist'); hist(dists(dists < 10^3));
    votes{i, 1} = find(dists < thresh*thresh);
    if (length(votes{i, 1}) > 2)
        %disp('hi');
    end
    if (length(votes{i, 1}) > length(votes{max_ind, 1}))
        max_ind = i;
        if (length(votes{i, 1}) > 0.9*num_pairs)
            break;
        end
    end
end
inlier_ind = votes{max_ind, 1};

tps1 = est_tps(pts1(inlier_ind, :), pts2(inlier_ind, 1));
tps2 = est_tps(pts1(inlier_ind, :), pts2(inlier_ind, 2));
tps = [tps1;tps2];

end