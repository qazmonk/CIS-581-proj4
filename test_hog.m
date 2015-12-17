%%
VISUALS = false;
%%SETUP REFERENCE FACE
grid_size = 8;
face = imread('poses_scaled/Nate/1202151736b.jpg');
mask = imread('poses_scaled/Nate/1202151736b-mask.png');
video = 'easy2';
[face_nr, face_nc, ~] = size(face);
frame_names = dir(sprintf('%s/frame*.png', video));
[face_xs, face_ys] = meshgrid(1:grid_size:face_nc, 1:grid_size: ...
                              face_nr);
face_inds = sub2ind([face_nr, face_nc], face_ys(:), face_xs(:));
face_inds(mask(face_inds) == 0) = [];
[face_ys, face_xs] = ind2sub([face_nr, face_nc], face_inds);

[face_features, face_points] = extractFeatures(face, ...
                                               [face_xs, face_ys]);
%This is how we 'disable' weights. We just set them to all ones so
%basically no weighting
feature_dim = size(face_features, 2);
face_weights = ones(length(face_features), feature_dim);

%calculate center of face and shift origin of all face points to it
s = regionprops(mask, 'centroid');
center = s.Centroid;
face_vecs = face_points - repmat(center, length(face_points), 1);

%%FIND FACE IN EACH FRAME
for fidx = 1:length(frame_names)
    frame_name = sprintf('%s/frame-%03d.png', video, fidx);
    frame = imread(frame_name);

    [frame_nr, frame_nc, ~] = size(frame);


    [frame_xs, frame_ys] = meshgrid(1:grid_size:frame_nc, 1:grid_size:frame_nr);



    %%GET FEATURES AND GENERATE SCORE MAP
    [frame_features, frame_points] = extractFeatures(frame,...
                                                     [frame_xs(:), frame_ys(:)]);


    [ score_map, voters, IDS, D ] = generate_score_map(frame_features, frame_points,...
                                                      face_features, face_vecs,...
                                                      face_weights,...
                                                      frame_nr, frame_nc);
    
    if VISUALS
        figure('Name', 'Score Map'); imagesc(score_map); axis image;
    end

    
    %%FIND MOST VOTED FOR POINT
    [~, max_inds] = max(score_map(:));
    [max_y, max_x] = ind2sub([frame_nr frame_nc], max_inds);
    votes = reshape(voters(max_y,  max_x, :), ...
                    length(frame_points), 1);

    %hold over from when features could cast multiple votes
    did_vote = any(votes, 2);
    
    %Do some Matlabing to get all the appropriate features from the
    %did_vote vector
    matched_frame_points = frame_points(did_vote, :);
    frame_matched_idxs = find(did_vote);
    matched_weights = D(frame_matched_idxs);
    face_matched_idxs = IDS(frame_matched_idxs);
    matched_face_points = face_points(face_matched_idxs, :);

    num_matches = length(matched_frame_points);
    %shift all points so that the origin is the predicted center
    %we do this so that the origin means the same thing in both the
    %face and the frame space
    shifted_mframe_points = matched_frame_points - repmat([max_x,max_y], ...
                                                      num_matches, ...
                                                      1);
    shifted_mface_points = matched_face_points - repmat(center, ...
                                                      num_matches, ...
                                                      1);
    
    %RANSAC IT!
    if (num_matches < 5)
        fprintf('only found %d matches\n', num_matches);
    end
    [tps, inlier_ids] = ransac_est_tps(shifted_mframe_points, ...
                                       shifted_mface_points, ...
                                       4);
    
    %More Matlabing to get the inliers
    inlier_weights = matched_weights(inlier_ids);
    inlier_frame_points = matched_frame_points(inlier_ids, :);
    inlier_face_points = matched_face_points(inlier_ids, :);
    [inlier_weights, sorted_inds] = sort(inlier_weights);
    inlier_frame_points = inlier_frame_points(sorted_inds, :);
    inlier_face_points = inlier_face_points(sorted_inds, :);

    %%ENFORE ONE-2-ONE MATCHING BEFORE WARPING
    %when multiple feature points match one face point we retain
    %the match with the lowest weight
    [o2o_face_points, o2o_inds, ~] = unique(inlier_face_points, ...
                                            'rows', 'first');
    o2o_frame_points = inlier_frame_points(o2o_inds, :);

    %Do the replacement
    replaced = replace_face(frame, face, mask, [max_x max_y], ...
                            center, o2o_frame_points, ...
                            o2o_face_points);
    output_name = sprintf('output2/frame-%03d.png', fidx);
    fprintf('writing to %s\n', output_name);
    imwrite(replaced, output_name);
   
    
    if VISUALS
        figure(); imagesc(replaced); axis image;
        figure(); imagesc(frame); axis image; hold on;
        Xin = matched_frame_points(inlier_ids, 1);
        Yin = matched_frame_points(inlier_ids, 2);
        X = matched_frame_points(:, 1); 
        Y = matched_frame_points(:, 2);
        U = face_vecs(face_matched_idxs(inlier_ids), 1);
        V = face_vecs(face_matched_idxs(inlier_ids), 2);
        for n = 1:length(X)
            th = 0:pi/50:2*pi;
            xpts = grid_size*cos(th)+X(n);
            ypts = grid_size*sin(th)+Y(n);
            if (any(inlier_ids == n))
                plot(xpts, ypts, 'Color', 'green');
            else
                plot(xpts, ypts, 'Color', 'blue');
            end

        end
        quiver(Xin, Yin, -U, -V, 0);
        hold off;
    end
                %
                %saveas(h, output_name);

end
