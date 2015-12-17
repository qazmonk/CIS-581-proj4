function [IDS, D] = weighted_knnsearch( wpts, ws, pts)

num_wpts = size(wpts, 1);
num_pts = size(pts, 1);
dim = size(ws, 2);

%make a version of each pts point that has been dotted with each of
%the weights (i.e dotted_pts(i, j, :) contains frame point i dotted
%with weight vector j
dotted_pts = zeros(num_pts, num_wpts, dim);
for i = 1:num_wpts
    dotted_pts(:, i, :) = repmat(ws(i, :), num_pts, 1) .* pts;
end

%compute distances and find the minimum
D = zeros(num_pts, 1);
IDS = zeros(num_pts, 1);
for i = 1:num_pts
    shaped_dotted_pts = reshape(dotted_pts(i, :, :), num_wpts, ...
                                dim);
    dists = sqrt(sum((shaped_dotted_pts - wpts) .^ 2, 2));
    [D(i), IDS(i)] = min(dists);
end

end