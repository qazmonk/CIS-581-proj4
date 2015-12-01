function [ y, x, rmax ] = anms( cimg, max_pts)
%ANMS Summary of this function goes here
%   Detailed explanation goes here
[counts, edges] = hist(cimg(:), 0:0.001:0.1);
thresh = 0.0005;
%thresh = edges(find(cumsum(counts) > numel(cimg) - 5000, 1))+0.0005;
[ y, x ] = find(cimg > thresh);
num_points = length(y);
while (num_points < max_pts)
    thresh = thresh*0.8;
    [ y, x ] = find(cimg > thresh);
    num_points = length(y);
end
dists = pdist2([x, y], [x, y]).^2;
values = cimg(sub2ind(size(cimg), y, x));
nearest_greater = [zeros(num_points, 1), y, x];
for i = 1:num_points
    tmp = min(dists(i, values > values(i)));
    if (~isempty(tmp))
        nearest_greater(i, 1) = tmp;
    else
        nearest_greater(i, 1) = Inf;
    end
end
rmax = max(nearest_greater(:, 1));
nearest_greater(nearest_greater == Inf) = rmax;
nearest_greater = flipud(sortrows(nearest_greater));

rmin_idx = find(nearest_greater(:, 1) == nearest_greater(1), 1, 'last');
next_rmin_idx = find(nearest_greater(:, 1) ==...
                     nearest_greater(rmin_idx+1, 1), 1, 'last');

while next_rmin_idx < max_pts
    rmin_idx = next_rmin_idx;
    next_rmin_idx = find(nearest_greater(:, 1) == ...
                         nearest_greater(rmin_idx+1, 1), 1, 'last');
end
rmax= sqrt(nearest_greater(rmin_idx, 1));
y = nearest_greater(1:rmin_idx, 2);
x = nearest_greater(1:rmin_idx, 3);
end
