%%

grid_size = 8;
face = imread('poses_scaled/Nate/1202151736b.jpg');
mask = imread('poses_scaled/Nate/1202151736b-mask.png');
[face_nr, face_nc, ~] = size(face);
frame_names = dir('easy3/frame*.png');
[face_xs, face_ys] = meshgrid(1:grid_size:face_nc, 1:grid_size: ...
                              face_nr);
face_inds = sub2ind([face_nr, face_nc], face_ys(:), face_xs(:));
face_inds(mask(face_inds) == 0) = [];
[face_ys, face_xs] = ind2sub([face_nr, face_nc], face_inds);

[face_features, face_points] = extractFeatures(face, ...
                                               [face_xs, face_ys]);
s = regionprops(mask, 'centroid');
center = s.Centroid;
face_vecs = face_points - repmat(center, length(face_points), 1);
for fidx = 1:length(frame_names)
    fprintf('computing frame %03d\n', fidx);
    frame_name = sprintf('easy1/frame-%03d.png', fidx);
    frame = imread(frame_name);

    [frame_nr, frame_nc, ~] = size(frame);


    [frame_xs, frame_ys] = meshgrid(1:grid_size:frame_nc, 1:grid_size:frame_nr);




    [frame_features, frame_points] = extractFeatures(frame,...
                                                     [frame_xs(:), frame_ys(:)]);



    
    %figure('Name', 'Features'); imagesc(face); axis image; hold on;
    %plot(face_vis, 'Color', 'green');

    %%

    K = 1;
    [IDS, D] = knnsearch(face_features, frame_features, 'K', K);
    score_map = zeros(frame_nr, frame_nc);
    [Xs, Ys] = meshgrid(1:frame_nc, 1:frame_nr);
    perc = 0;
    prev_perc = 0;

    voters = false(frame_nr, frame_nc, length(frame_features), K);
    for i = 1:length(frame_features)

        for k = 1:K
            w = exp(-D(i, k));

            if (w > 0.3)
                pred_cent = frame_points(i, :) - face_vecs(IDS(i, k), :);
                if (pred_cent(1) > 0 && pred_cent(1) < frame_nc &&...
                    pred_cent(2) > 0 && pred_cent(2) < frame_nr)
                    
                    xmin = max(1,  floor(pred_cent(1))-12);
                    xmax = min(frame_nc, floor(pred_cent(1))+12);
                    ymin = max(1, floor(pred_cent(2))-12);
                    ymax = min(frame_nr, floor(pred_cent(2))+12);

                    
                    ydist = w*normpdf(ymin:ymax, pred_cent(2), 8)';
                    xdist = w*normpdf(xmin:xmax, pred_cent(1), 8);
                    
                    voters(ymin:ymax, xmin:xmax, i, k) = true;
                    dist = conv2(ydist, xdist, 'full');
                    score_map(ymin:ymax, xmin:xmax) = (score_map(ymin:ymax, ...
                                                                 xmin:xmax) ...
                                                       + dist);
                end
                %dist = reshape(mvnpdf([Ys(:), Xs(:)], pred_cent, [1, 1]), [frame_nr, frame_nc]);

                %score_map = score_map + dist;    
            end
        end
    end % for i = 1:length(frame_features)

    [~, max_inds] = max(score_map(:));
    [max_y, max_x] = ind2sub([frame_nr frame_nc], max_inds);
    h=figure('Name', 'Voters'); set(gcf,'Visible', 'off'); imagesc(frame); axis image; hold on;

    for m = 1:length(max_inds) 
        votes = reshape(voters(max_y(m),  max_x(m), :, :), ...
                        [length(frame_points) K]);

        
        did_vote = any(votes, 2);
        
        matched_frame_points = frame_points(did_vote, :);
        frame_matched_idxs = find(did_vote);
        matched_weights = D(frame_matched_idxs);
        face_matched_idxs = IDS(frame_matched_idxs);
        matched_face_points = face_points(face_matched_idxs, :);

        num_matches = length(matched_frame_points);
        shifted_mframe_points = matched_frame_points - repmat([max_x,max_y], ...
                                                          num_matches, ...
                                                          1);
        shifted_mface_points = matched_face_points - repmat(center, ...
                                                          num_matches, ...
                                                          1);
        [tps, inlier_ids] = ransac_est_tps(shifted_mframe_points, ...
                                           shifted_mface_points, ...
                                           4);
        
        inlier_weights = matched_weights(inlier_ids);
        Xin = matched_frame_points(inlier_ids, 1);
        Yin = matched_frame_points(inlier_ids, 2);
        X = matched_frame_points(:, 1);
        Y = matched_frame_points(:, 2);
        U = face_vecs(matches(inlier_ids), 1);
        V = face_vecs(matches(inlier_ids), 2);
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
        output_name = sprintf('output2/frame-%03d.png', fidx);
        saveas(h, output_name);
    end
end